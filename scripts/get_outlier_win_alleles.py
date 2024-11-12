# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3
import tempfile
import os


#create d stats table
def create_windows_table(conn, win_size, pop, outlier_type, win_allele_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {win_allele_table}")
  cur.execute(f"""create table {win_allele_table}
                  (odwsa_id int primary key,
                   win_id varchar(30),
                   chrom varchar(20),
                   pos int,
                   sample_id varchar(50),
                   num_der_alleles int)""")

                   
  cur.execute(f"create index idx_odwsa_pos_{win_size}{pop}{outlier_type} on {win_allele_table}(pos)")
  cur.execute(f"create index idx_odwsa_dsw_{win_size}{pop}{outlier_type} on {win_allele_table}(win_id)")
  cur.execute(f"create index idx_odwsa_sample_id_{win_size}{pop}{outlier_type} on {win_allele_table}(sample_id)")

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




def main():
  win_size = int(sys.argv[1])
  pop = sys.argv[2]
  outlier_type = str(sys.argv[3])
  pop_str = str(sys.argv[4])
  zarr_prefix = str(sys.argv[5])
  db_file = str(sys.argv[6])
  
  
  
  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  d_win_table = "d_stat_win_" + str(win_size) + '_' + pop_str
  poly_win_table = "poly_win_" + str(win_size)
  sites_table = "outlier_" + outlier_type + "_win_sites_" + str(win_size) + '_' + pop_str
  win_allele_table = "outlier_" + outlier_type + "_win_alleles_" + pop + "_" + str(win_size) + '_' + pop_str
  
  create_windows_table(conn, win_size, pop, outlier_type, win_allele_table)
  
  
  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_species s, sample_pop_link l
                            where s.sample_id = l.sample_id""", conn)

  outlier_sql = f"""select win_id, chrom, start, end
                    from {d_win_table}
                    where win_id in (select distinct win_id from {sites_table})"""
                    
  outlier_win_df = pd.read_sql(outlier_sql, conn)
  
  
  #get outgroup
  outgroup = pop_str.split('_')[3]
  
  
  # Intialize pop dictionary.
  idx_pop_dicc = {}
  for pop_i in list(set(meta_df['pop'])):
    idx_pop_dicc[pop_i] = meta_df['vcf_order'][meta_df['pop'] == pop_i]

  
  # Get alleles for the user provided pop for every outlier window 
  sample_ids = list(meta_df['sample_id'][meta_df['pop'] == pop])
  odwsa_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t:
    for idx in range(outlier_win_df.shape[0]):
      win_id = outlier_win_df.win_id.values[idx]
      chrom = outlier_win_df.chrom.values[idx]
      start = outlier_win_df.start.values[idx]
      end = outlier_win_df.end.values[idx]
      print(win_id, chrom, start, end)
       
      # Extract the genotype callset and positions.
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos = load_callset_pos(chrom, zarr_file)
      # Identify the window to extract.
      wind_loc = all_pos.locate_range(start, end)
      win_gt = allel.GenotypeArray(callset[wind_loc])
      
      
      pop_gt=win_gt.take(idx_pop_dicc[pop], axis=1)
      outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1)
      
      total_alleles = pop_gt.count_alleles().sum(axis=1)
      # If there are no altenative alleles...
      if (pop_gt.count_alleles().shape[1] == 1):
          # Calculate alternative allele frequencies.
          alt_alleles = pop_gt.count_alleles()[:, 0] - 1
      else:
          # Calculate alternative allele frequencies.
          alt_alleles = pop_gt.count_alleles()[:, 1]    
                
      for pos_i in range(pop_gt.shape[0]):
        if alt_alleles[pos_i] not in [total_alleles[pos_i], 0]:
          pos = all_pos[wind_loc][pos_i]
          
          for sample_i, sample_id in enumerate(sample_ids):
            num_alt_alleles = pop_gt[pos_i, sample_i][0] + pop_gt[pos_i, sample_i][1]  
            num_der_alleles = num_alt_alleles
            if outgroup_gt[pos_i,][0][0] == 1 and num_alt_alleles != -2:
              num_der_alleles = 2 - num_alt_alleles
            odwsa_id += 1
            t.write(f"{odwsa_id}\t{win_id}\t{chrom}\t{pos}\t{sample_id}\t{num_der_alleles}\n")  
        
          
    t.flush()         
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} {win_allele_table}" """)

           
  

if __name__ == '__main__':
  main()
