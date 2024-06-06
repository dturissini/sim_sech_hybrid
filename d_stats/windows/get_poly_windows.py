# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create d stats table
def create_windows_table(conn, win_size, poly_win_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {poly_win_table}")
  cur.execute(f"""create table {poly_win_table}
                  (pw_id int primary key,
                   chrom varchar(20),
                   start int,
                   end int,
                   num_sites int,
                   seg_sites int,
                   der_freq_sim float,
                   der_freq_sech float,
                   der_freq_ssh float,
                   der_freq_diff_sim_sech float,
                   der_freq_diff_sim_ssh float,
                   der_freq_diff_sech_ssh float,
                   pi_sim float,
                   pi_sech float,
                   pi_ssh float,
                   dxy_sim_sech float,
                   dxy_sim_ssh float,
                   dxy_sech_ssh float,
                   fst_sim_sech float,
                   fst_sim_ssh float,
                   fst_sech_ssh float)""")

  cur.execute(f"create index idx_pw_start_end_{win_size} on {poly_win_table}(start,end)")
  cur.execute(f"create index idx_pw_end_{win_size} on {poly_win_table}(end)")

  cur.execute(f"drop table if exists ssh_der_freq")
  cur.execute(f"""create table ssh_der_freq
                  (sdf_id int primary key,
                   pop varchar(4),
                   der_freq decimal(3,2),
                   num_sites int)""")
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


def calc_mean_der_freq(gt, outgroup_gt, sfs):
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs = gt.count_alleles().to_frequencies()[:, 1]    

    der_freqs = alt_freqs 
    for i, freq in enumerate(der_freqs):  
      if outgroup_gt[i,][0][0] == 1:
        der_freqs[i] = 1 - freq  
 
      rounded_der_freq = round(der_freqs[i], 2)
      if not np.isnan(rounded_der_freq):
        if rounded_der_freq in sfs:
          sfs[rounded_der_freq] += 1
        else:
          sfs[rounded_der_freq] = 1
    
    idx_nz = [i for i in range(len(der_freqs)) if der_freqs[i] > 0 and der_freqs[i] < 1] 
    mean_der_freq = np.nansum(der_freqs[idx_nz]) / der_freqs[idx_nz].size
    
    return mean_der_freq, sfs


def calc_mean_der_freq_diff(gt_a, gt_b, outgroup_gt):
    # If there are no altenative alleles...
    if (gt_a.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs_a = gt_a.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs_a = gt_a.count_alleles().to_frequencies()[:, 1]    

    der_freqs_a = alt_freqs_a
    for i, freq in enumerate(der_freqs_a):   
      if outgroup_gt[i,][0][0] == 1:
        der_freqs_a[i] = 1 - freq    

    # If there are no altenative alleles...
    if (gt_b.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_freqs_b = gt_b.count_alleles().to_frequencies()[:, 0] - 1
    # Else...
    else:
        # Calculate alternative allele frequencies.
        alt_freqs_b = gt_b.count_alleles().to_frequencies()[:, 1]    

    der_freqs_b = alt_freqs_b
    for i, freq in enumerate(der_freqs_b):   
      if outgroup_gt[i,][0][0] == 1:
        der_freqs_b[i] = 1 - freq    
    
    idx_nz = [i for i in range(len(der_freqs_a)) if der_freqs_a[i] > 0 and der_freqs_a[i] < 1 and der_freqs_b[i] > 0 and der_freqs_b[i] < 1] 
    
    der_freq_diff = abs(der_freqs_a[idx_nz] - der_freqs_b[idx_nz])
    mean_der_freq_diff = np.nansum(der_freq_diff) / der_freq_diff[~np.isnan(der_freq_diff)].size
    
    return mean_der_freq_diff


# Define a function to calculate nucleotide diversity.
def pixy(gt, pop_idx):
    # Extract the alternative allele count matrix.
    aac = gt.take(pop_idx, axis=1).count_alleles()
    # If the alternative allele count matrix is mono-allelic.
    if aac.shape[1] == 1:
        # Then there are no pairwise differences.
        return 0
    # Else, the site is bi-allelic.
    else:
        # Determine the number of called chromosomes for each site.
        called_chromosomes = np.nansum(aac, axis=1)
        # Create a mask where there are no called genotypes.
        mask = called_chromosomes == 0
        # Determine the allele counts of the derived/alternative allele.
        derived_allele_count = aac[:, 1]
        # Determine the allele counts of the ancestral/reference allele.
        ancestral_allele_count = called_chromosomes - derived_allele_count
        # Determine the number of comparisons per site.
        nC2 = np.array([((n * (n - 1)) / 2) for n in called_chromosomes])
        # Calculate the numerator.
        numerator = np.nansum((derived_allele_count[~mask] * ancestral_allele_count[~mask]))
        # Calculate the denominator.
        denominator = np.nansum(nC2[~mask])
        # Calculate pixy.
        pi_pixy = numerator / denominator
        return pi_pixy


# Define a function to calculate the dXY for a given locus.
def calc_dxy(gt, pop_x, pop_y):
    # Compute the allele frequencies.
    pop_x_freqs = calc_alt_freqs(gt=gt.take(pop_x, axis=1))
    pop_y_freqs = calc_alt_freqs(gt=gt.take(pop_y, axis=1))
    # Calculate the per site dXY.
    per_site_dxy = ((pop_x_freqs * (1 - pop_y_freqs)) + (pop_y_freqs * (1 - pop_x_freqs)))
    # Calculate the average dXY for this locus.
    dxy = np.nansum(per_site_dxy) / per_site_dxy[~np.isnan(per_site_dxy)].size
    return dxy


# Define a function to calculate fst.
def calc_fst(gt, pop_a, pop_b):
    # Determine allele counts.
    a_ac = gt.take(pop_a, axis=1).count_alleles()
    b_ac = gt.take(pop_b, axis=1).count_alleles()
    # Calculate the numerator and denominator for Hudson's Fst estimator.
    a_b_num, a_b_den = allel.hudson_fst(a_ac, b_ac)
    # Calculate Fst.
    a_b_fst = np.nansum(a_b_num) / np.nansum(a_b_den)
    return a_b_fst


def main():
  win_size = int(sys.argv[1])
  zarr_prefix = str(sys.argv[2])
  db_file = str(sys.argv[3])
  
  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  poly_win_table = "poly_win_" + str(win_size)
  create_windows_table(conn, win_size, poly_win_table)  

  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql('select sample, pop from sample_pop', conn)
  
  # Intialize pop dictionary.
  idx_dicc = {}
  for pop in ['simulans', 'sim_sech_hybrid', 'sechellia', 'melanogaster']:
      # Fill the dictionary.
      idx_dicc[pop] = meta_df[meta_df['pop'] == pop].index.values

  
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


  pw_ids = []
  der_freq_sims = []
  der_freq_sechs = []
  der_freq_sshs = []
  der_freq_diff_sim_sechs = []
  der_freq_diff_sim_sshs = []
  der_freq_diff_sech_sshs = []
  pi_sims = []
  pi_sechs = []
  pi_sshs = []
  dxy_sim_sechs = []
  dxy_sim_sshs = []
  dxy_sech_sshs = []
  fst_sim_sechs = []
  fst_sim_sshs = []
  fst_sech_sshs = [] 
  
  sfs_sim = {}
  sfs_sech = {}
  sfs_ssh = {}
  # For every ortholog window.
  for idx in range(window_df.shape[0]):
      pw_ids.append(idx + 1)
      # Extract the ortholog information.
      chrom = chroms[idx]
      start = starts[idx]
      end = ends[idx]
      num_sites = tot_sites[idx]
      num_seg_sites = seg_sites[idx]
       
      if num_sites == 0:
        der_freq_sims.append(np.nan)
        der_freq_sechs.append(np.nan)
        der_freq_sshs.append(np.nan)
        der_freq_diff_sim_sechs.append(np.nan)
        der_freq_diff_sim_sshs.append(np.nan)
        der_freq_diff_sech_sshs.append(np.nan)
        pi_sims.append(np.nan)
        pi_sechs.append(np.nan)
        pi_sshs.append(np.nan)
        dxy_sim_sshs.append(np.nan)
        dxy_sim_sechs.append(np.nan)
        dxy_sech_sshs.append(np.nan)
        fst_sim_sechs.append(np.nan)
        fst_sim_sshs.append(np.nan)
        fst_sech_sshs.append(np.nan)
      else:
        # Extract the genotype callset and positions.
        zarr_file = zarr_prefix + '_' + chrom + '.zarr'
        callset, all_pos = load_callset_pos(chrom, zarr_file)
        # Identify the window to extract.
        wind_loc = all_pos.locate_range(start, end)
        
        win_gt = allel.GenotypeArray(callset[wind_loc])
        
        # Compute derived allele freq
        der_freq_sim, sfs_sim = calc_mean_der_freq(gt=win_gt.take(idx_dicc['simulans'], axis=1), outgroup_gt=win_gt.take(idx_dicc['melanogaster'], axis=1), sfs=sfs_sim)
        der_freq_sech, sfs_sech = calc_mean_der_freq(gt=win_gt.take(idx_dicc['sechellia'], axis=1), outgroup_gt=win_gt.take(idx_dicc['melanogaster'], axis=1), sfs=sfs_sech)
        der_freq_ssh, sfs_ssh = calc_mean_der_freq(gt=win_gt.take(idx_dicc['sim_sech_hybrid'], axis=1), outgroup_gt=win_gt.take(idx_dicc['melanogaster'], axis=1), sfs=sfs_ssh)
        
        # Compute derived allele freq diff        
        der_freq_diff_sim_sech = calc_mean_der_freq_diff(gt_a=win_gt.take(idx_dicc['simulans'], axis=1), gt_b=win_gt.take(idx_dicc['sechellia'], axis=1), outgroup_gt=win_gt.take(idx_dicc['melanogaster'], axis=1))
        der_freq_diff_sim_ssh = calc_mean_der_freq_diff(gt_a=win_gt.take(idx_dicc['simulans'], axis=1), gt_b=win_gt.take(idx_dicc['sim_sech_hybrid'], axis=1), outgroup_gt=win_gt.take(idx_dicc['melanogaster'], axis=1))
        der_freq_diff_sech_ssh = calc_mean_der_freq_diff(gt_a=win_gt.take(idx_dicc['sechellia'], axis=1), gt_b=win_gt.take(idx_dicc['sim_sech_hybrid'], axis=1), outgroup_gt=win_gt.take(idx_dicc['melanogaster'], axis=1))
        
        # Compute pi
        pi_sim = pixy(gt=win_gt, pop_idx=idx_dicc['simulans'])
        pi_sech = pixy(gt=win_gt, pop_idx=idx_dicc['sechellia'])
        pi_ssh = pixy(gt=win_gt, pop_idx=idx_dicc['sim_sech_hybrid'])

        # Compute dxy
        dxy_sim_sech = calc_dxy(gt=win_gt, pop_x=idx_dicc['simulans'], pop_y=idx_dicc['sechellia'])
        dxy_sim_ssh = calc_dxy(gt=win_gt, pop_x=idx_dicc['simulans'], pop_y=idx_dicc['sim_sech_hybrid'])
        dxy_sech_ssh = calc_dxy(gt=win_gt, pop_x=idx_dicc['sechellia'], pop_y=idx_dicc['sim_sech_hybrid'])

        # Compute Fst
        fst_sim_sech = calc_fst(gt=win_gt, pop_a=idx_dicc['simulans'], pop_b=idx_dicc['sechellia'])
        fst_sim_ssh = calc_fst(gt=win_gt, pop_a=idx_dicc['simulans'], pop_b=idx_dicc['sim_sech_hybrid'])
        fst_sech_ssh = calc_fst(gt=win_gt, pop_a=idx_dicc['sechellia'], pop_b=idx_dicc['sim_sech_hybrid'])


        der_freq_sims.append(der_freq_sim)
        der_freq_sechs.append(der_freq_sech)
        der_freq_sshs.append(der_freq_ssh)
        der_freq_diff_sim_sechs.append(der_freq_diff_sim_sech)
        der_freq_diff_sim_sshs.append(der_freq_diff_sim_ssh)
        der_freq_diff_sech_sshs.append(der_freq_diff_sech_ssh)
        pi_sims.append(pi_sim)
        pi_sechs.append(pi_sech)
        pi_sshs.append(pi_ssh)
        dxy_sim_sechs.append(dxy_sim_sech)
        dxy_sim_sshs.append(dxy_sim_ssh)
        dxy_sech_sshs.append(dxy_sech_ssh)                                                    
        fst_sim_sechs.append(fst_sim_sech)
        fst_sim_sshs.append(fst_sim_ssh)
        fst_sech_sshs.append(fst_sech_ssh)
           
  
  #add columns to window_df
  window_df.insert(0, 'pw_id', pw_ids)
  
  window_df["der_freq_sim"] = der_freq_sims
  window_df["der_freq_sech"] = der_freq_sechs
  window_df["der_freq_ssh"] = der_freq_sshs
  window_df["der_freq_diff_sim_sech"] = der_freq_diff_sim_sechs
  window_df["der_freq_diff_sim_ssh"] = der_freq_diff_sim_sshs
  window_df["der_freq_diff_sech_ssh"] = der_freq_diff_sech_sshs
  window_df["pi_sim"] = pi_sims
  window_df["pi_sech"] = pi_sechs
  window_df["pi_ssh"] = pi_sshs
  window_df["dxy_sim_sech"] = dxy_sim_sechs
  window_df["dxy_sim_ssh"] = dxy_sim_sshs
  window_df["dxy_sech_ssh"] = dxy_sech_sshs
  window_df["fst_sim_sech"] = fst_sim_sechs
  window_df["fst_sim_ssh"] = fst_sim_sshs
  window_df["fst_sech_ssh"] = fst_sech_sshs
  
  
  #import window_df to db
  conn = sqlite3.connect(db_file)
  window_df.to_sql(poly_win_table, if_exists = 'append', index=False, con=conn)

  sdf_id = 0
  for freq in sfs_sim:
    sdf_id += 1
    conn.execute(f"""insert into ssh_der_freq
                     values
                     ({sdf_id}, 'sim', {freq}, {sfs_sim[freq]})""")
    
  for freq in sfs_sech:
    sdf_id += 1
    conn.execute(f"""insert into ssh_der_freq
                     values
                     ({sdf_id}, 'sech', {freq}, {sfs_sech[freq]})""")
    
  for freq in sfs_ssh:
    sdf_id += 1
    conn.execute(f"""insert into ssh_der_freq
                     values
                     ({sdf_id}, 'ssh', {freq}, {sfs_ssh[freq]})""")
    
  conn.commit()
  conn.close()
  

if __name__ == '__main__':
  main()
