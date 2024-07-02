# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import os
import tempfile
import sqlite3


# sys.argv[0] = reference
# sys.argv[1] = win_size
# sys.argv[2] = zarr_prefix
# sys.argv[3] = windows_file


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

# Extract the ref assembly.
ref = str(sys.argv[1])
win_size = int(sys.argv[2])
zarr_prefix = str(sys.argv[3])
db_file = str(sys.argv[4])

# Intialize a dictionary of chromosome lengths.
if ref == 'sim':
    chromosome_dicc = {
        '2L': 23857595, 
        '2R': 22319025,
        '3L': 23399903, 
        '3R': 28149585,
        'X': 22032822,
    }
    
    
# Compute the adjusted chromosome lengths for windowing.
adj_chrom_dicc = chr_seq_len(win_size, chromosome_dicc)

# Intialize a dictionary to store the results.
df_dicc = {
    'chr': [],
    'start': [],
    'end': [],
    'tot_sites': [],
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
            df_dicc['chr'].append(chrom)
            df_dicc['start'].append(left)
            df_dicc['end'].append(right)
            df_dicc['tot_sites'].append(wind_idx.size)
            df_dicc['seg_sites'].append(var_mask.sum())
        # Else, there are no sites in this window.
        else:
            # Fill the dictionary.
            df_dicc['chr'].append(chrom)
            df_dicc['start'].append(left)
            df_dicc['end'].append(right)
            df_dicc['tot_sites'].append(0)
            df_dicc['seg_sites'].append(0)


conn = sqlite3.connect(db_file)  
cur = conn.cursor()    
                        
cur.execute("drop table if exists windows")
cur.execute("""create table windows
               (w_id int primary key,
                chrom varchar(20),
                start int,
                end int,
                num_sites int,
                seg_sites int)""")
 
cur.execute("create index idx_w_start_end on windows(start,end)")
cur.execute("create index idx_w_end on windows(end)")
conn.close()

with tempfile.NamedTemporaryFile(mode='w') as o: 
  w_id = 0
  for i, chrom in enumerate(df_dicc['chr']):
    w_id += 1
    o.write(f"{w_id}\t'{chrom}'\t{df_dicc['start'][i]}\t{df_dicc['end'][i]}\t{df_dicc['tot_sites'][i]}\t{df_dicc['seg_sites'][i]}\n")


  o.flush() 
  os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {o.name} windows" """)



### Convert the dictionary to a dataframe.
##df = pd.DataFrame(df_dicc)
##df_filtered = df[df['tot_sites'] > 0]
##
### Export the window information.
##df_filtered.to_csv(win_file, sep='\t', index=False)