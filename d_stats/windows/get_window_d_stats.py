# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create d stats table
def create_windows_table(conn, win_size, table_suffix, d_stats_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {d_stats_table}")
  cur.execute(f"""create table {d_stats_table}
                  (win_id varchar(30) primary key,
                   chrom varchar(20),
                   start int,
                   end int,
                   num_sites int,
                   seg_sites int,
                   abba float,
                   baba float,
                   baaa float,
                   abaa float,
                   d float,
                   d_anc float,
                   d_plus float)""")

  cur.execute(f"create index idx_w_start_end_{win_size}_{table_suffix} on {d_stats_table}(start,end)")
  cur.execute(f"create index idx_w_end_{win_size}_{table_suffix} on {d_stats_table}(end)")
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
        start_stop[index] = [window_start, (window_start+window_size)-1] # Sci-kit allele uses [inclusive, inclsuive] indexing.
        index += 1
    return windows, start_stop



#calculate alternative allele frequencies
def calc_alt_freqs(gt):
    if (gt.count_alleles().shape[1] == 1):
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    else:
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs


#get abba-baba site patterns
def site_patterns(p1, p2, p3, p4):
    # Calculate site pattern counts.
    abba = np.nansum((1 - p1) * (p2) * (p3) * (1 - p4))
    baba = np.nansum((p1) * (1 - p2) * (p3) * (1 - p4))
    baaa = np.nansum((p1) * (1 - p2) * (1 - p3) * (1 - p4))
    abaa = np.nansum((1 - p1) * (p2) * (1 - p3) * (1 - p4))
    return abba, baba, baaa, abaa


#calculate site patterns
def dros_site_patterns(
    gt,
    p1_idx, p2_idx, p3_idx, p4_idx,
):
    #determine the indicies where each pop has called genotypes
    p1_mask = (gt.take(p1_idx, axis=1).is_called() == True).any(axis=1)
    p2_mask = (gt.take(p2_idx, axis=1).is_called() == True).any(axis=1)
    p3_mask = (gt.take(p3_idx, axis=1).is_called() == True).any(axis=1)
    p4_mask = (gt.take(p4_idx, axis=1).is_called() == True).any(axis=1)

    #get the indices where all pops have called genotypes
    called_mask = (p1_mask & p2_mask & p3_mask & p4_mask)

    if (called_mask.sum() == 0):
        abba, baba, baaa, abaa = np.zeros(4)
    else:
        #get the indicies with polymorphic sites
        var_mask = gt.compress(called_mask, axis=0).count_alleles().is_variant()
        
        if (var_mask.sum() == 0):
            abba, baba, baaa, abaa = np.zeros(4)
        # Else...
        else:
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p4_alt_freqs = calc_alt_freqs(gt.take(p4_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))

            #polarize the allele frequencies based on the most common allele in the outgroup
            p1_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            p4_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p4_alt_freqs - 1), p4_alt_freqs)

            #calculate the site pattern counts.
            abba, baba, baaa, abaa = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs, p4_der_freqs,
            )
    return abba, baba, baaa, abaa

#get d stats
def get_d_stats(abba, baba, baaa, abaa):
  d_num = abba - baba
  d_den = abba + baba
  d_anc_num = baaa - abaa
  d_anc_den = baaa + abaa
  d_plus_num = d_num + d_anc_num
  d_plus_den = d_den + d_anc_den
  
  d = np.nan
  d_anc = np.nan
  d_plus = np.nan
  
  if d_den != 0:
    d = d_num / float(d_den)
  
  if d_anc_den != 0:
    d_anc = d_anc_num / float(d_anc_den)
  
  if d_plus_den != 0:
    d_plus = d_plus_num / float(d_plus_den)
  
  return d, d_anc, d_plus


def main():
  p1 = str(sys.argv[1])
  p2 = str(sys.argv[2])
  p3 = str(sys.argv[3])
  p4 = str(sys.argv[4])
  win_size = int(sys.argv[5])
  zarr_prefix = str(sys.argv[6])
  db_file = str(sys.argv[7])
  
  table_suffix = p1 + '_' + p2 + '_' + p3 + '_' + p4
    
    
  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  d_stats_table = "d_stat_win_" + str(win_size) + '_' + table_suffix
  create_windows_table(conn, win_size, table_suffix, d_stats_table)  
  
  
  #make pandas dataframe of metadata
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_species s, sample_pop_link l
                            where s.sample_id = l.sample_id
                            and l.pop in ('{p1}', '{p2}', '{p3}', '{p4}')""", conn)
  
  #intialize the pop dictionary
  idx_pop_dicc = {}
  for pop in list(set(meta_df['pop'])):
    idx_pop_dicc[pop] = meta_df['vcf_order'][meta_df['pop'] == pop]

  
  #intialize a dictionary of chromosome lengths
  chrom_query = conn.execute("""select chrom, chrom_len 
                                from chrom_lens""")
  chromosome_dicc = {}
  for chrom, chrom_len in chrom_query:
    chromosome_dicc[chrom] = chrom_len
      
      
  #get adjusted chromosome lengths for windowing
  adj_chrom_dicc = chr_seq_len(win_size, chromosome_dicc)
  
  #intialize a dictionary to store the window results
  df_dicc = {
      'win_id': [],
      'chrom': [],
      'start': [],
      'end': [],
      'num_sites': [],
      'seg_sites': [],
  }

  #make dictionary of windows
  for chrom in adj_chrom_dicc:
      #extract the genotype callset and positions
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
  
    
  #convert the dictionary to a dataframe
  window_df = pd.DataFrame(df_dicc)  
  
  
  #get d stat values
  win_ids = window_df.win_id.values
  chroms = window_df.chrom.values
  starts = window_df.start.values
  ends = window_df.end.values
  tot_sites = window_df.num_sites.values
  seg_sites = window_df.seg_sites.values
  
  abbas = []
  babas = []
  abaas = []
  baaas = []
  ds = []
  d_ancs = []
  d_pluses = [] 
  for idx in range(window_df.shape[0]):
      chrom = chroms[idx]
      start = starts[idx]
      end = ends[idx]
      num_sites = tot_sites[idx]
      num_seg_sites = seg_sites[idx]

       
      if num_sites == 0:
        abbas.append(np.nan)
        babas.append(np.nan)
        baaas.append(np.nan)
        abaas.append(np.nan)
        ds.append(np.nan)
        d_ancs.append(np.nan)
        d_pluses.append(np.nan)
      else:
        zarr_file = zarr_prefix + '_' + chrom + '.zarr'
        callset, all_pos = load_callset_pos(chrom, zarr_file)
        wind_loc = all_pos.locate_range(start, end)
        abba, baba, baaa, abaa = dros_site_patterns(gt=allel.GenotypeArray(callset[wind_loc]),
                                                    p1_idx=idx_pop_dicc[p1], p2_idx=idx_pop_dicc[p2],
                                                    p3_idx=idx_pop_dicc[p3], p4_idx=idx_pop_dicc[p4],
                                                    )
        abbas.append(abba)
        babas.append(baba)
        baaas.append(baaa)
        abaas.append(abaa)
        
        d, d_anc, d_plus = get_d_stats(abba, baba, baaa, abaa)                                                        
        ds.append(d)
        d_ancs.append(d_anc)
        d_pluses.append(d_plus)
           
  
  #add lists to window_df
  window_df["abba"] = abbas
  window_df["baba"] = babas
  window_df["baaa"] = baaas
  window_df["abaa"] = abaas
  window_df["d"] = ds
  window_df["d_anc"] = d_ancs
  window_df["d_plus"] = d_pluses
  
  
  #import window_df to db
  conn = sqlite3.connect(db_file)  
  window_df.to_sql(d_stats_table, if_exists = 'append', index=False, con=conn)
  conn.close()

if __name__ == '__main__':
  main()
