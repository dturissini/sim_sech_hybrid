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
def create_tables(inv_conn, sech_allele_table):
  cur = inv_conn.cursor()    
  cur.execute(f"drop table if exists {sech_allele_table}")
  cur.execute(f"""create table {sech_allele_table}
                  (isa_id int primary key,
                   inv_name varchar(20),
                   chrom varchar(20),
                   pos int,
                   sample_id varchar(50),
                   num_der_alleles int,
                   num_anro_alleles int,
                   outgroup_allele varchar(1),
                   der_allele varchar(1),
                   anro_allele varchar(1))""")

                   
  cur.execute(f"create index idx_isa_pos on {sech_allele_table}(pos)")
  cur.execute(f"create index idx_isa_sample_id on {sech_allele_table}(sample_id)")

  cur.close()
  


# Define a function to load genotyope and positions arrays.
def load_callset_pos(chrom, zarr_file):
    # Load the vcf file.
    callset = zarr.open_group(zarr_file, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    refs = callset[f'{chrom}/variants/REF']
    alts = callset[f'{chrom}/variants/ALT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos, refs, alts




def main():
  zarr_prefix = str(sys.argv[1])
  d_win_db_file = str(sys.argv[2])
  inv_db_file = str(sys.argv[3])
  inversion_file = str(sys.argv[4])
  anro_line = str(sys.argv[5])
  
  
  
  
  #establish db connection and create d_stats_win table
  d_win_conn = sqlite3.connect(d_win_db_file)  
  inv_conn = sqlite3.connect(inv_db_file)  
  
  
  sech_allele_table = "inversion_sech_alleles"
  
  create_tables(inv_conn, sech_allele_table)
  
  
  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_species s, sample_pop_link l
                            where s.sample_id = l.sample_id""", d_win_conn)

  #process and store inversion locations
  inversions = {}
  with open(inversion_file, 'r') as i:
    next(i)
    for line in i:
      inv_name, chrom, start, end = line.strip().split()
      inversions[inv_name] = {'chrom': chrom, 'start': start, 'end': end}


  
  #define pops
  pop = 'sech'
  outgroup = 'mel'
  
  
  # Intialize pop dictionary.
  idx_pop_dicc = {}
  for pop_i in list(set(meta_df['pop'])):
    idx_pop_dicc[pop_i] = meta_df['vcf_order'][meta_df['pop'] == pop_i]
  
  anro_idx = meta_df['vcf_order'][meta_df['sample_id'] == anro_line]
  
  # Get alleles for sechellia
  sample_ids = list(meta_df['sample_id'][meta_df['pop'] == pop])
  isa_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t:
    for inv_name in inversions:
      chrom = inversions[inv_name]['chrom']
      start = int(inversions[inv_name]['start'])
      end = int(inversions[inv_name]['end'])
      print(inv_name, chrom, start, end)
       
      # Extract the genotype callset and positions.
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos, refs, alts = load_callset_pos(chrom, zarr_file)
      # Identify the window to extract.
      wind_loc = all_pos.locate_range(start, end)
      win_gt = allel.GenotypeArray(callset[wind_loc])
      
      refs_win = refs[wind_loc]
      alts_win = alts[wind_loc]
      
      
      pop_gt=win_gt.take(idx_pop_dicc[pop], axis=1)
      outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1)
      anro_gt=win_gt.take(anro_idx, axis=1)
      
      outgroup_alleles_numeric = outgroup_gt[:, 0, 0]
      
      total_alleles = pop_gt.count_alleles().sum(axis=1)
      # If there are no alternative alleles...
      if (pop_gt.count_alleles().shape[1] == 1):
          # Calculate alternative allele frequencies.
          alt_alleles = pop_gt.count_alleles()[:, 0] - 1
      else:
          # Calculate alternative allele frequencies.
          alt_alleles = pop_gt.count_alleles()[:, 1]    
      
      outgroup_allele_i = 0
      for pos_i in range(pop_gt.shape[0]):    
        alleles = [refs_win[pos_i][0], alts_win[pos_i][0]]   
        der_allele_numeric = abs(outgroup_alleles_numeric[outgroup_allele_i] - 1)
        anro_allele_numeric = anro_gt[pos_i,][0][0]
        
        outgroup_allele = alleles[outgroup_alleles_numeric[outgroup_allele_i]]
        der_allele = alleles[der_allele_numeric]
        anro_allele = alleles[anro_allele_numeric]
        outgroup_allele_i += 1
        
        if alt_alleles[pos_i] not in [total_alleles[pos_i], 0]:                    
          if anro_gt[pos_i,][0][0] != -1:
            pos = all_pos[wind_loc][pos_i]
                        
            for sample_i, sample_id in enumerate(sample_ids):
              num_alt_alleles = pop_gt[pos_i, sample_i][0] + pop_gt[pos_i, sample_i][1]  
              
              if num_alt_alleles >= 0:
                num_der_alleles = num_alt_alleles
                if outgroup_gt[pos_i,][0][0] == 1 and num_alt_alleles != -2:
                  num_der_alleles = 2 - num_alt_alleles
                
                num_anro_alleles = 0
                if anro_gt[pos_i,][0][0] == 0:
                  num_anro_alleles = 2 - num_alt_alleles
                
                isa_id += 1
                t.write(f"{isa_id}\t{inv_name}\t{chrom}\t{pos}\t{sample_id}\t{num_der_alleles}\t{num_anro_alleles}\t{outgroup_allele}\t{der_allele}\t{anro_allele}\n")  
        
          
    t.flush()         
    os.system(f"""sqlite3 {inv_db_file} ".mode tabs" ".import {t.name} {sech_allele_table}" """)

           
  

if __name__ == '__main__':
  main()
