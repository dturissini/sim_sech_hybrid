# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3



def create_windows_table(conn, win_size, poly_win_table, sfs_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {poly_win_table}")
  cur.execute(f"""create table {poly_win_table}
                  (pw_id int primary key,
                   win_id varchar(30),
                   chrom varchar(20),
                   start int,
                   end int,
                   num_sites int,
                   seg_sites int,
                   pop varchar(50),
                   der_freq float,
                   pi float)""")
  
  cur.execute(f"create index idx_pw_start_end_{win_size} on {poly_win_table}(start,end)")
  cur.execute(f"create index idx_pw_end_{win_size} on {poly_win_table}(end)")
  cur.execute(f"create index idx_pw_pop_{win_size} on {poly_win_table}(pop)")
  cur.execute(f"create index idx_pw_win_id_{win_size} on {poly_win_table}(win_id)")
    
  
  cur.execute(f"drop table if exists {sfs_table}")
  cur.execute(f"""create table {sfs_table}
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


def calc_mean_der_freq(gt, outgroup_gt, sfs_win):
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
        if rounded_der_freq in sfs_win:
          sfs_win[rounded_der_freq] += 1
        else:
          sfs_win[rounded_der_freq] = 1
    
    idx_nz = [i for i in range(len(der_freqs)) if der_freqs[i] > 0 and der_freqs[i] < 1] 
    if der_freqs[idx_nz].size > 0:
      mean_der_freq = np.nansum(der_freqs[idx_nz]) / der_freqs[idx_nz].size
    else:
      mean_der_freq = np.nan
    
    return mean_der_freq, sfs_win



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



def main():
  win_size = int(sys.argv[1])
  outgroup = str(sys.argv[2])
  zarr_prefix = str(sys.argv[3])
  db_file = str(sys.argv[4])
  


  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  poly_win_table = "poly_win_" + str(win_size)
  sfs_table = "sfs_der_freq_" + str(win_size)
  
  create_windows_table(conn, win_size, poly_win_table, sfs_table)

  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_pop s, lk_pop l
                            where s.sample_id = l.sample_id""", conn)
  
  pops = list(set(meta_df['pop']))
  
  # Intialize pop dictionary.
  idx_dicc = {}
  for pop in pops:
    idx_dicc[pop] = meta_df['vcf_order'][meta_df['pop'] == pop]
  
  #remove outgroups from pop since we don't want to calculate polymorphism measures for it
  pops.remove(outgroup)
  
  # Intialize a dictionary of chromosome lengths.
  chrom_query = conn.execute("""select chrom, chrom_len 
                                from chrom_lens""")
  chromosome_dicc = {}
  for chrom, chrom_len in chrom_query:
    chromosome_dicc[chrom] = chrom_len
      
      
  # Compute the adjusted chromosome lengths for windowing.
  adj_chrom_dicc = chr_seq_len(win_size, chromosome_dicc)
  
  # Intialize a dictionary to store the results.
             
  # For every chromosome.
  pw_id = 0
  sdf_id = 0
  sfs= {}
  for pop in pops:
    sfs[pop] = {}
    
  for chrom in adj_chrom_dicc:
      # Extract the genotype callset and positions.
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos = load_callset_pos(chrom, zarr_file)
      # Construct the window dictionaries.
      wind_dicc, left_right = window_info(
          all_pos, win_size, adj_chrom_dicc[chrom],
      )
      
      for pop in pops:
        print(chrom, pop)
        df_dicc = {'win_id': [],
                   'chrom': [],
                   'start': [],
                   'end': [],
                   'num_sites': [],
                   'seg_sites': []}
        
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
                df_dicc['win_id'].append(chrom + '_' + str(left))
                df_dicc['chrom'].append(chrom)
                df_dicc['start'].append(left)
                df_dicc['end'].append(right)
                df_dicc['num_sites'].append(wind_idx.size)
                df_dicc['seg_sites'].append(var_mask.sum())
            # Else, there are no sites in this window.
            else:
                # Fill the dictionary.
                df_dicc['win_id'].append(chrom + '_' + str(left))
                df_dicc['chrom'].append(chrom)
                df_dicc['start'].append(left)
                df_dicc['end'].append(right)
                df_dicc['num_sites'].append(0)
                df_dicc['seg_sites'].append(0)
  
    
        # Convert the dictionary to a dataframe.
        window_df = pd.DataFrame(df_dicc)  
        
        # Extract ortholog information.
        win_ids = window_df.win_id.values
        chroms = window_df.chrom.values
        starts = window_df.start.values
        ends = window_df.end.values
        tot_sites = window_df.num_sites.values


        pw_ids = []
        table_pops = []
        der_freqs = []
        pis = []
  
        # For every ortholog window.
        for idx in range(window_df.shape[0]):
            pw_id += 1
            pw_ids.append(pw_id)
            table_pops.append(pop)
            # Extract the ortholog information.
            chrom = chroms[idx]
            start = starts[idx]
            end = ends[idx]
            num_sites = tot_sites[idx]
            
            if num_sites == 0:
              der_freqs.append(np.nan)
              pis.append(np.nan)
            else:
              # Extract the genotype callset and positions.
              zarr_file = zarr_prefix + '_' + chrom + '.zarr'
              callset, all_pos = load_callset_pos(chrom, zarr_file)
              # Identify the window to extract.
              wind_loc = all_pos.locate_range(start, end)
              
              win_gt = allel.GenotypeArray(callset[wind_loc])
              
              # Compute derived allele freq
              der_freq, sfs[pop] = calc_mean_der_freq(gt=win_gt.take(idx_dicc[pop], axis=1), outgroup_gt=win_gt.take(idx_dicc[outgroup], axis=1), sfs_win=sfs[pop])
                      
              # Compute pi
              pi = pixy(gt=win_gt, pop_idx=idx_dicc[pop])
        
        
        
              der_freqs.append(der_freq)
              pis.append(pi)
                 
        
        #add columns to window_df
        window_df.insert(0, 'pw_id', pw_ids)
        
        window_df["pop"] = table_pops
        window_df["der_freq"] = der_freqs
        window_df["pi"] = pis
        
        
        #import window_df to db
        conn = sqlite3.connect(db_file)
        window_df.to_sql(poly_win_table, if_exists = 'append', index=False, con=conn)
 
 
  for pop in sfs:      
    for freq in sfs[pop]:
      sdf_id += 1
      conn.execute(f"""insert into {sfs_table}
                       values
                       ({sdf_id}, '{pop}', {freq}, {sfs[pop][freq]})""")
        
  conn.commit()
  conn.close()
  

if __name__ == '__main__':
  main()
