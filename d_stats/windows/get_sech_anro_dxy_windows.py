# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create database table
def create_windows_table(conn, win_size, dxy_win_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {dxy_win_table}")
  cur.execute(f"""create table {dxy_win_table}
                  (sadw_id int primary key,
                   win_id varchar(30),
                   chrom varchar(20),
                   start int,
                   end int,
                   num_sites int,
                   seg_sites int,
                   sample_id varchar(50),
                   dxy_sech_anro float)""")

  cur.execute(f"create index idx_sadw_start_end_{win_size} on {dxy_win_table}(start,end)")
  cur.execute(f"create index idx_sadw_end_{win_size} on {dxy_win_table}(end)")
  cur.execute(f"create index idx_sadw_sample_id_{win_size} on {dxy_win_table}(sample_id)")
  cur.execute(f"create index idx_sadw_win_id_{win_size} on {dxy_win_table}(win_id)")
  cur.close()
  


#load genotyope and positions arrays
def load_callset_pos(chrom, zarr_file):
    callset = zarr.open_group(zarr_file, mode='r')
    geno = callset[f'{chrom}/calldata/GT']
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos

#compute adjusted chromosome lengths
def chr_seq_len(window_size, chr_dicc):
    new_chr_dicc = {}
    for key in chr_dicc :
        chr_len = chr_dicc[key]
        new_chr_len = (chr_len//window_size)*window_size
        new_chr_dicc[key] = new_chr_len
    return new_chr_dicc

#break up a chromosome into windows
def window_info(positions, window_size, sequence_length):
    windows = {}
    start_stop = {}
    index = 0
    for window_start in range(1, int(sequence_length), int(window_size)):
        windows[index] = np.where(((window_start <= positions) & (positions < (window_start+window_size))))[0]
        start_stop[index] = [window_start, (window_start+window_size)-1] # Sci-kit allele uses [inclusive, inclsuive] indexing
        index += 1
    return windows, start_stop


#calculate alternative allele frequencies
def calc_alt_freqs(gt):
    if (gt.count_alleles().shape[1] == 1):
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    else:
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs


#calculate dxy for a window
def calc_dxy(gt, pop_x, pop_y):
    sparse_x = gt.take(pop_x, axis=1).to_sparse()
    sparse_y = gt.take(pop_y, axis=1).to_sparse()
    
    if -1 * len(sparse_x.data) == sum(sparse_x.data) or -1 * len(sparse_y.data) == sum(sparse_y.data):
      dxy=np.nan
    else:
      pop_x_freqs = calc_alt_freqs(gt=gt.take(pop_x, axis=1))
      pop_y_freqs = calc_alt_freqs(gt=gt.take(pop_y, axis=1))
      per_site_dxy = ((pop_x_freqs * (1 - pop_y_freqs)) + (pop_y_freqs * (1 - pop_x_freqs)))
      dxy = np.nansum(per_site_dxy) / per_site_dxy[~np.isnan(per_site_dxy)].size
    return dxy



def main():
  win_size = int(sys.argv[1])
  zarr_prefix = str(sys.argv[2])
  db_file = str(sys.argv[3])
      

  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  dxy_win_table = "sech_anro_dxy_win_" + str(win_size)
  create_windows_table(conn, win_size, dxy_win_table)  

  #make a pandas dataframe with metadata
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_species s, sample_pop_link l
                            where s.sample_id = l.sample_id""", conn)
  
  #idenfity the indices for the samples of interest
  focal_sample_ids = meta_df['sample_id'][meta_df['pop'] in ['sech', 'sechden', 'sechlad', 'sechmari', 'sechpras', 'ssh']]
  
  #intialize dictionary with pop indices
  idx_pop_dicc = {}
  for pop in list(set(meta_df['pop'])):
    idx_dicc[pop] = meta_df['vcf_order'][meta_df['pop'] == pop]
    
      

  
  #make a dictionary of chromosome lengths
  chrom_query = conn.execute("""select chrom, chrom_len 
                                from chrom_lens""")
  chromosome_dicc = {}
  for chrom, chrom_len in chrom_query:
    chromosome_dicc[chrom] = chrom_len
      
      
  #get the adjusted chromosome lengths for windowing
  adj_chrom_dicc = chr_seq_len(win_size, chromosome_dicc)
  
  #intialize a dictionary to store the results
  df_dicc = {
      'win_id': [],
      'chrom': [],
      'start': [],
      'end': [],
      'num_sites': [],
      'seg_sites': [],
  }
  
  # get stats for each window
  for chrom in adj_chrom_dicc:
      #get the genotype callset and positions
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos = load_callset_pos(chrom, zarr_file)

      wind_dicc, left_right = window_info(
          all_pos, win_size, adj_chrom_dicc[chrom],
      )
      for wind in wind_dicc:
          left, right = left_right[wind]
          wind_idx = np.where(((left <= all_pos) & (all_pos <= right)))[0]

          if wind_idx.size > 0:
              wind_loc = all_pos.locate_range(left, right)
              sub_gt = allel.GenotypeArray(callset[wind_loc])
              var_mask = sub_gt.count_alleles().is_variant()

              df_dicc['win_id'].append(chrom + '_' + str(left))
              df_dicc['chrom'].append(chrom)
              df_dicc['start'].append(left)
              df_dicc['end'].append(right)
              df_dicc['num_sites'].append(wind_idx.size)
              df_dicc['seg_sites'].append(var_mask.sum())
          else:
              df_dicc['win_id'].append(chrom + '_' + str(left))
              df_dicc['chrom'].append(chrom)
              df_dicc['start'].append(left)
              df_dicc['end'].append(right)
              df_dicc['num_sites'].append(0)
              df_dicc['seg_sites'].append(0)
  
    
  #get dxy for each sample
  sadw_id = 0
  for sample_id in focal_sample_ids:
    print(sample_id)

    window_df = pd.DataFrame(df_dicc)  
    
    # Extract ortholog information.
    win_ids = window_df.win_id.values
    chroms = window_df.chrom.values
    starts = window_df.start.values
    ends = window_df.end.values
    tot_sites = window_df.num_sites.values
    seg_sites = window_df.seg_sites.values
    
    
    sadw_ids = []
    sample_ids = []
    dxy_sech_anros = []
    
    for idx in range(window_df.shape[0]):
        sadw_id += 1
        sadw_ids.append(sadw_id)
        sample_ids.append(sample_id)
        chrom = chroms[idx]
        start = starts[idx]
        end = ends[idx]
        num_sites = tot_sites[idx]
        num_seg_sites = seg_sites[idx]
         
        if num_sites == 0:
          dxy_sech_anros.append(np.nan)
        else:
          zarr_file = zarr_prefix + '_' + chrom + '.zarr'
          callset, all_pos = load_callset_pos(chrom, zarr_file)

          wind_loc = all_pos.locate_range(start, end)
          
          win_gt = allel.GenotypeArray(callset[wind_loc])
                  
    
          #compute dxy
          dxy_sech_anro = calc_dxy(gt=win_gt, pop_x=idx_dicc['sech_anro'], pop_y=idx_dicc[sample_id])
    
          dxy_sech_anros.append(dxy_sech_anro)                                                 
             
    
    #add columns to window_df
    window_df.insert(0, 'sadw_id', sadw_ids)  
    window_df["sample_id"] = sample_ids
    window_df["dxy_sech_anro"] = dxy_sech_anros
    
    
    #import window_df to db
    conn = sqlite3.connect(db_file)
    window_df.to_sql(dxy_win_table, if_exists = 'append', index=False, con=conn)
    
      
    conn.commit()
    conn.close()
  

if __name__ == '__main__':
  main()
