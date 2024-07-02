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
                  (dsw_id int primary key,
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
  


# Define a function to load genotyope and positions arrays.
def load_callset_pos(chrom, zarr_file):
    # Load the vcf file.
    callset = zarr.open_group(zarr_file, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos

# Define a function to compute adjusted chromosome lengths.
def chr_seq_len(window_size, chr_dicc):
    # Initialize new empty dictionary to store the output.
    new_chr_dicc = {}
    # Iterate through every chromosome.
    for key in chr_dicc :
        # Floor divide every chromosome length by the window size and multiply
        # to find the new chromosome length.
        chr_len = chr_dicc[key]
        new_chr_len = (chr_len//window_size)*window_size
        # Refill dictionary with the new chromosome lengths.
        new_chr_dicc[key] = new_chr_len
    return new_chr_dicc

# Define a function to break up a chromosome into windows.
def window_info(positions, window_size, sequence_length):
    # Intialize a dicctionary with the start and stop position for each window.
    windows = {}
    start_stop = {}
    index = 0
    for window_start in range(1, int(sequence_length), int(window_size)):
        windows[index] = np.where(((window_start <= positions) & (positions < (window_start+window_size))))[0]
        start_stop[index] = [window_start, (window_start+window_size)-1] # Sci-kit allele uses [inclusive, inclsuive] indexing.
        index += 1
    return windows, start_stop



# Define a function to calculate alternative allele frequencies.
def calc_alt_freqs(gt):
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]
    return alt_freqs

# Define a site pattern function.
def site_patterns(p1, p2, p3, p4):
    # Calculate site pattern counts.
    abba = np.nansum((1 - p1) * (p2) * (p3) * (1 - p4))
    baba = np.nansum((p1) * (1 - p2) * (p3) * (1 - p4))
    baaa = np.nansum((p1) * (1 - p2) * (1 - p3) * (1 - p4))
    abaa = np.nansum((1 - p1) * (p2) * (1 - p3) * (1 - p4))
    return abba, baba, baaa, abaa

# Define a function to calculate site patterns.
def dros_site_patterns(
    gt,
    p1_idx, p2_idx, p3_idx, p4_idx,
):
    # Determine the indicies where each population has called genotypes.
    p1_mask = (gt.take(p1_idx, axis=1).is_called() == True).any(axis=1)
    p2_mask = (gt.take(p2_idx, axis=1).is_called() == True).any(axis=1)
    p3_mask = (gt.take(p3_idx, axis=1).is_called() == True).any(axis=1)
    p4_mask = (gt.take(p4_idx, axis=1).is_called() == True).any(axis=1)
    # Determine the indicied where all populations have called genotypes.
    called_mask = (p1_mask & p2_mask & p3_mask & p4_mask)
    # If there are no sites called between all three samples...
    if (called_mask.sum() == 0):
        # Set the results to np.nan since we don't have any sites to perform computations on.
        abba, baba, baaa, abaa = np.zeros(6)
    # Else...
    else:
        # Determine the indicies where we have varibale sites.
        var_mask = gt.compress(called_mask, axis=0).count_alleles().is_variant()
        # If there are no variable sites...
        if (var_mask.sum() == 0):
            # Set the results to 0 since we are iterating over QC'ed regions.
            abba, baba, baaa, abaa = np.zeros(6)
        # Else...
        else:
            # Calculate the alternative allele frequencies.
            p1_alt_freqs = calc_alt_freqs(gt.take(p1_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p2_alt_freqs = calc_alt_freqs(gt.take(p2_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p3_alt_freqs = calc_alt_freqs(gt.take(p3_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            p4_alt_freqs = calc_alt_freqs(gt.take(p4_idx, axis=1).compress(called_mask, axis=0).compress(var_mask, axis=0))
            # Polarize the allele frequencies based on the most common allele in the outgroup.
            p1_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p1_alt_freqs - 1), p1_alt_freqs)
            p2_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p2_alt_freqs - 1), p2_alt_freqs)
            p3_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p3_alt_freqs - 1), p3_alt_freqs)
            p4_der_freqs = np.where(p4_alt_freqs > 0.5, np.abs(p4_alt_freqs - 1), p4_alt_freqs)
            # Calculate the site pattern counts.
            abba, baba, baaa, abaa = site_patterns(
                p1_der_freqs, p2_der_freqs, p3_der_freqs, p4_der_freqs,
            )
    return abba, baba, baaa, abaa


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
  
  
  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql(f"""select s.sample_id, subpop 
                            from sample_pop s, lk_subpop l
                            where s.sample_id = l.sample_id
                            and subpop in ('{p1}', '{p2}', '{p3}', '{p4}')""", conn)
  
  # Intialize pop dictionary.
  idx_dicc = {}
  for pop in [p1, p2, p3, p4]:
      # Fill the dictionary.
      idx_dicc[pop] = meta_df[meta_df['subpop'] == pop].index.values

  
  # Intialize a dictionary of chromosome lengths.
  chrom_query = conn.execute("""select chrom, chrom_len 
                                from chrom_lens""")
  chromosome_dicc = {}
  for chrom, chrom_len in chrom_query:
    chromosome_dicc[chrom] = chrom_len
      
      
  # Compute the adjusted chromosome lengths for windowing.
  adj_chrom_dicc = chr_seq_len(win_size, chromosome_dicc)
  
  # Intialize a dictionary to store the results.
  df_dicc = {
      'chrom': [],
      'start': [],
      'end': [],
      'num_sites': [],
      'seg_sites': [],
  }
  # For every chromosome.
  for chrom in adj_chrom_dicc:
      # Extract the genotype callset and positions.
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos = load_callset_pos(chrom, zarr_file)
      # Construct the window dictionaries.
      wind_dicc, left_right = window_info(
          all_pos, win_size, adj_chrom_dicc[chrom],
      )
      # For every window.
      for wind in wind_dicc:
          # Extract the left and right positions.
          left, right = left_right[wind]
          # Determine the what sites are in this region.
          wind_idx = np.where(((left <= all_pos) & (all_pos <= right)))[0]
          # If the window has at least one site.
          if wind_idx.size > 0:
              # Identify the window to extract.
              wind_loc = all_pos.locate_range(left, right)
              # Subset the genotype matrix.
              sub_gt = allel.GenotypeArray(callset[wind_loc])
              # Determine which positions are segregating.
              var_mask = sub_gt.count_alleles().is_variant()
              # Fill the dictionary.
              df_dicc['chrom'].append(chrom)
              df_dicc['start'].append(left)
              df_dicc['end'].append(right)
              df_dicc['num_sites'].append(wind_idx.size)
              df_dicc['seg_sites'].append(var_mask.sum())
          # Else, there are no sites in this window.
          else:
              # Fill the dictionary.
              df_dicc['chrom'].append(chrom)
              df_dicc['start'].append(left)
              df_dicc['end'].append(right)
              df_dicc['num_sites'].append(0)
              df_dicc['seg_sites'].append(0)
  
    
  # Convert the dictionary to a dataframe.
  window_df = pd.DataFrame(df_dicc)  
  
  
  
  # Extract ortholog information.
  chroms = window_df.chrom.values
  starts = window_df.start.values
  ends = window_df.end.values
  tot_sites = window_df.num_sites.values
  seg_sites = window_df.seg_sites.values
  
  dsw_ids = []
  abbas = []
  babas = []
  abaas = []
  baaas = []
  ds = []
  d_ancs = []
  d_pluses = [] 
  # For every ortholog window.
  for idx in range(window_df.shape[0]):
      dsw_ids.append(idx + 1)
      # Extract the ortholog information.
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
        # Extract the genotype callset and positions.
        zarr_file = zarr_prefix + '_' + chrom + '.zarr'
        callset, all_pos = load_callset_pos(chrom, zarr_file)
        # Identify the window to extract.
        wind_loc = all_pos.locate_range(start, end)
        # Calculate site patterns.
        abba, baba, baaa, abaa = dros_site_patterns(gt=allel.GenotypeArray(callset[wind_loc]),
                                                    p1_idx=idx_dicc[p1], p2_idx=idx_dicc[p2],
                                                    p3_idx=idx_dicc[p3], p4_idx=idx_dicc[p4],
                                                    )
        abbas.append(abba)
        babas.append(baba)
        baaas.append(baaa)
        abaas.append(abaa)
        
        d, d_anc, d_plus = get_d_stats(abba, baba, baaa, abaa)                                                        
        ds.append(d)
        d_ancs.append(d_anc)
        d_pluses.append(d_plus)
           
  
  #add columns to window_df
  window_df.insert(0, 'dsw_id', dsw_ids)
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
