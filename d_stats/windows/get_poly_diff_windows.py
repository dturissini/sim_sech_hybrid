# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create poly_diff_win_table
def create_windows_table(conn, win_size, poly_diff_win_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {poly_diff_win_table}")
  cur.execute(f"""create table {poly_diff_win_table}
                  (pwd_id int primary key,
                   win_id varchar(30),
                   chrom varchar(20),
                   start int,
                   end int,
                   num_sites int,
                   seg_sites int,
                   pop_a varchar(50),
                   pop_b varchar(50),
                   der_freq_diff float,
                   dxy float,
                   fst float)""")
  
  cur.execute(f"create index idx_pwd_start_end_{win_size} on {poly_diff_win_table}(start,end)")
  cur.execute(f"create index idx_pwd_end_{win_size} on {poly_diff_win_table}(end)")
  cur.execute(f"create index idx_pwd_pops_{win_size} on {poly_diff_win_table}(pop_a, pop_b)")
  cur.execute(f"create index idx_pwd_win_id_{win_size} on {poly_diff_win_table}(win_id)")
  cur.close()



#define a function to load genotyope and positions arrays
def load_callset_pos(chrom, zarr_file):
    # Load the vcf file.
    callset = zarr.open_group(zarr_file, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos


#define a function to compute adjusted chromosome lengths
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


#define a function to break a chromosome up into windows
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



#calculate alternative allele frequencies
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


#calculate the difference for derived allele frequencies between two pops
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
    if der_freq_diff[~np.isnan(der_freq_diff)].size > 0:
      mean_der_freq_diff = np.nansum(der_freq_diff) / der_freq_diff[~np.isnan(der_freq_diff)].size
    else:
      mean_der_freq_diff = np.nan
    
    return mean_der_freq_diff



#calculate dxy for a given window
def calc_dxy(gt, pop_x, pop_y):
    # Compute the allele frequencies.
    pop_x_freqs = calc_alt_freqs(gt=gt.take(pop_x, axis=1))
    pop_y_freqs = calc_alt_freqs(gt=gt.take(pop_y, axis=1))
    # Calculate the per site dXY.
    per_site_dxy = ((pop_x_freqs * (1 - pop_y_freqs)) + (pop_y_freqs * (1 - pop_x_freqs)))
    # Calculate the average dXY for this locus.
    dxy = np.nansum(per_site_dxy) / per_site_dxy[~np.isnan(per_site_dxy)].size
    return dxy


#calculate Fst for a given window
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
  outgroup = str(sys.argv[2])
  zarr_prefix = str(sys.argv[3])
  db_file = str(sys.argv[4])
  


  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  poly_diff_win_table = "poly_diff_win_" + str(win_size)
  
  create_windows_table(conn, win_size, poly_diff_win_table)

  #make a pandas dataframe of metadata
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_species s, sample_pop_link l
                            where s.sample_id = l.sample_id""", conn)
  
  pops = list(set(meta_df['species']))
  pops.sort()
  
  
  #intialize a pop dictionary
  idx_dicc = {}
  for pop in pops:
    idx_dicc[pop] = meta_df['vcf_order'][meta_df['species'] == pop]
  
  #remove outgroups from pop since we don't want to calculate polymorphism measures for it
  pops.remove(outgroup)
  
  #intialize a dictionary of chromosome lengths
  chrom_query = conn.execute("""select chrom, chrom_len 
                                from chrom_lens""")
  chromosome_dicc = {}
  for chrom, chrom_len in chrom_query:
    chromosome_dicc[chrom] = chrom_len
      
      
  #compute the adjusted chromosome lengths for windowing
  adj_chrom_dicc = chr_seq_len(win_size, chromosome_dicc)
  
  #create a dictionary to store the results             
  pwd_id = 0
  sdf_id = 0
  for chrom in adj_chrom_dicc:
      #get the genotype callset and positions
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos = load_callset_pos(chrom, zarr_file)
      #make the window dictionaries
      wind_dicc, left_right = window_info(
          all_pos, win_size, adj_chrom_dicc[chrom],
      )
      
      for j, pop_b in enumerate(pops):
        for i in range(j):
          pop_a = pops[i]
          #don't calculate differences between subpopulations
          if pop_a[:4] != pop_b[:4]:
            print(chrom, pop_a, pop_b)
            df_dicc = {'win_id': [],
                       'chrom': [],
                       'start': [],
                       'end': [],
                       'num_sites': [],
                       'seg_sites': []}
            
            for wind in wind_dicc:
                #extract the left and right boundaries of the window
                left, right = left_right[wind]
                #find the sites in this window
                wind_idx = np.where(((left <= all_pos) & (all_pos <= right)))[0]
                #process the window
                if wind_idx.size > 0:
                    #get the window location
                    wind_loc = all_pos.locate_range(left, right)
                    #subset the genotype matrix
                    sub_gt = allel.GenotypeArray(callset[wind_loc])
                    #find polymorphic sites
                    var_mask = sub_gt.count_alleles().is_variant()
                    #add to lists
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
            
            #make lists of window values
            win_ids = window_df.win_id.values
            chroms = window_df.chrom.values
            starts = window_df.start.values
            ends = window_df.end.values
            tot_sites = window_df.num_sites.values
            
            #make empty lists corresponding to database fields
            pwd_ids = []
            pop_as = []
            pop_bs = []
            der_freq_diffs = []
            dxys = []
            fsts = []
            
            #get stats for each window
            for idx in range(window_df.shape[0]):
                pwd_id += 1
                pwd_ids.append(pwd_id)
                pop_as.append(pop_a)
                pop_bs.append(pop_b)
                #extract the ortholog information
                chrom = chroms[idx]
                start = starts[idx]
                end = ends[idx]
                num_sites = tot_sites[idx]
                
                
                if num_sites == 0:
                  der_freq_diffs.append(np.nan)
                  dxys.append(np.nan)
                  fsts.append(np.nan)
                else:
                  #extract the genotypes and positions
                  zarr_file = zarr_prefix + '_' + chrom + '.zarr'
                  callset, all_pos = load_callset_pos(chrom, zarr_file)
                  #get the windows
                  wind_loc = all_pos.locate_range(start, end)
                  
                  win_gt = allel.GenotypeArray(callset[wind_loc])
                  
                  #compute derived allele freq difference
                  der_freq_diff = calc_mean_der_freq_diff(gt_a=win_gt.take(idx_dicc[pop_a], axis=1), gt_b=win_gt.take(idx_dicc[pop_b], axis=1), outgroup_gt=win_gt.take(idx_dicc[outgroup], axis=1))
                          
                  #compute dxy
                  dxy = calc_dxy(gt=win_gt, pop_x=idx_dicc[pop_a], pop_y=idx_dicc[pop_b])
                  
                  #compute fst
                  fst = calc_fst(gt=win_gt, pop_a=idx_dicc[pop_a], pop_b=idx_dicc[pop_b])
            
                  der_freq_diffs.append(der_freq_diff)
                  dxys.append(dxy)
                  fsts.append(fst)
                     
            
            #add columns to window_df
            window_df.insert(0, 'pwd_id', pwd_ids)
            
            window_df["pop_a"] = pop_as
            window_df["pop_b"] = pop_bs
            window_df["der_freq_diff"] = der_freq_diffs
            window_df["dxy"] = dxys
            window_df["fst"] = fsts
            
            
            #import window_df to db
            conn = sqlite3.connect(db_file)
            window_df.to_sql(poly_diff_win_table, if_exists = 'append', index=False, con=conn)
                          
            conn.commit()


  conn.close()
  

if __name__ == '__main__':
  main()
