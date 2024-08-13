# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create gwas allele freq table
def create_freq_table(conn):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists gwas_allele_freq")
  cur.execute(f"""create table gwas_allele_freq
                  (gaf_id int primary key,
                   chrom varchar(20),
                   pos int,
                   p float,
                   sim_freq decimal(5,4),
                   sech_freq decimal(5,4),
                   ssh_freq decimal(5,4))""")
  
  cur.execute(f"create index idx_gaf_pos on gwas_allele_freq(pos)")
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


def main():
  zarr_prefix = '/work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel'                                                                                                     
  db_file = '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
  
  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  
  create_freq_table(conn)
  
  conn_gwas = sqlite3.connect('/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db')
  
  
  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql('select sample_id, pop from sample_pop', conn)
  
  # Intialize pop dictionary.
  idx_dicc = {}
  for pop in ['simulans', 'sechellia', 'sim_sech_hybrid']:
      # Fill the dictionary.
      idx_dicc[pop] = meta_df[meta_df['pop'] == pop].index.values
  
  
  # Intialize a dictionary to store the results.
  freq_dicc = {
      'chrom': [],
      'pos': [],
      'p': [],
      'sim_freq': [],
      'sech_freq': [],
      'ssh_freq': [],
  }
  
  chrom_query = conn_gwas.execute("""select distinct chrom
                                    from gwas_snp
                                    where p < 1e-7
                                    and chrom not in ('4')""")
  
  for (chrom, ) in chrom_query:    
    print(chrom)                             
    gwas_df = pd.read_sql(f"""select pos, p
                              from gwas_snp
                              where chrom = '{chrom}'
                              and p < 1e-7
                              order by pos""", conn_gwas)
       
    # Extract the genotype callset and positions.
    zarr_file = zarr_prefix + '_' + chrom + '.zarr'
    callset, all_pos = load_callset_pos(chrom, zarr_file)
      
    gwas_pos, gwas_idx, all_pos_idx = np.intersect1d(gwas_df['pos'], all_pos, assume_unique=True, return_indices=True)
    
    # Calculate the alternative allele frequencies.
    gt=allel.GenotypeArray(callset[all_pos_idx])
    sim_alt_freqs = calc_alt_freqs(gt.take(idx_dicc['simulans'], axis=1))
    sech_alt_freqs = calc_alt_freqs(gt.take(idx_dicc['sechellia'], axis=1))
    ssh_alt_freqs = calc_alt_freqs(gt.take(idx_dicc['sim_sech_hybrid'], axis=1))
      
    # Fill the dictionary.
    freq_dicc['chrom'].extend([chrom] * len(gwas_pos))
    freq_dicc['pos'].extend(gwas_df['pos'][gwas_idx])
    freq_dicc['p'].extend(gwas_df['p'][gwas_idx])
    freq_dicc['sim_freq'].extend(sim_alt_freqs)
    freq_dicc['sech_freq'].extend(sech_alt_freqs)
    freq_dicc['ssh_freq'].extend(ssh_alt_freqs)
    
    
    
  # Convert the dictionary to a dataframe.
  freq_df = pd.DataFrame(freq_dicc)  
  gaf_ids = [x + 1 for x in list(range(len(freq_dicc['chrom'])))]
  freq_df.insert(0, 'gaf_id', gaf_ids)
  
  
  #import window_df to db
  conn = sqlite3.connect(db_file)  
  freq_df.to_sql('gwas_allele_freq', if_exists = 'append', index=False, con=conn)
  conn.close()

if __name__ == '__main__':
  main()
