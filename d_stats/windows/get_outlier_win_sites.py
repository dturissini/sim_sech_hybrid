# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create d stats table
def create_windows_table(conn, win_size, outlier_type, sites_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {sites_table}")
  cur.execute(f"""create table {sites_table}
                  (odwps_id int primary key,
                   win_id varchar(30),
                   chrom varchar(20),
                   pos int,
                   mel_allele varchar(1),
                   total_alleles_sim int,
                   total_alleles_ssh int,
                   total_alleles_sech int,
                   total_alleles_sech_anro int,
                   total_alleles_sech_denis int,
                   total_alleles_sech_ladigue int,
                   total_alleles_sech_marianne int,
                   total_alleles_sech_praslin int,
                   total_alleles_sech_unknown int,
                   der_alleles_sim int,
                   der_alleles_ssh int,
                   der_alleles_sech int,
                   der_alleles_sech_anro int,
                   der_alleles_sech_denis int,
                   der_alleles_sech_ladigue int,
                   der_alleles_sech_marianne int,
                   der_alleles_sech_praslin int,
                   der_alleles_sech_unknown int)""")

                   
  cur.execute(f"create index idx_odwps_pos_{win_size}{outlier_type} on {sites_table}(pos)")
  cur.execute(f"create index idx_odwps_dsw_{win_size}{outlier_type} on {sites_table}(win_id)")

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


# get derived allele counts for a given window
def get_der_allele_counts(gt, outgroup_gt):
    total_alleles = gt.count_alleles().sum(axis=1)
    # If there are no altenative alleles...
    if (gt.count_alleles().shape[1] == 1):
        # Calculate alternative allele frequencies.
        alt_alleles = gt.count_alleles()[:, 0] - 1
    else:
        # Calculate alternative allele frequencies.
        alt_alleles = gt.count_alleles()[:, 1]    

    der_alleles = alt_alleles 
    for i, num_der in enumerate(der_alleles):  
      if outgroup_gt[i,][0][0] == 1:
        der_alleles[i] = total_alleles[i] - num_der      
    
    return total_alleles, der_alleles


def main():
  win_size = int(sys.argv[1])
  outlier_type = str(sys.argv[2])
  pop_str = str(sys.argv[3])
  zarr_prefix = str(sys.argv[4])
  db_file = str(sys.argv[5])
  
  
  
  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  

  d_win_table = "d_stat_win_" + str(win_size) + '_' + pop_str
  poly_win_table = "poly_win_" + str(win_size)
  sites_table = "outlier_" + outlier_type + "_win_sites_" + str(win_size) + '_' + pop_str

  create_windows_table(conn, win_size, outlier_type, sites_table)  


  # Read in meta data as a pandas dataframe.
  meta_df = pd.read_sql(f"""select s.sample_id, l.pop, vcf_order
                            from sample_species s, sample_pop_link l
                            where s.sample_id = l.sample_id""", conn)
  
  #define sql to identify outlier windows depending on the type of outlier we're looking for
  if outlier_type == 'd_plus':
    d_plus = pd.read_sql(f"""select d_plus
                             from {d_win_table}
                             where num_sites > 1000""", conn)
      
    d_plus_quantiles = np.quantile(d_plus['d_plus'], [.001, .99])
    win_sql = f"""select win_id, chrom, start, end
                  from {d_win_table}
                  where num_sites > 1000
                  and (d_plus > {d_plus_quantiles[1]} or d_plus < {d_plus_quantiles[0]})"""
  elif outlier_type[:2] == 'pi':
    pop = outlier_type.split('_')[1]
    pi = pd.read_sql(f"""select pi
                         from {poly_win_table}
                         where num_sites > 1000
                         and pop = '{pop}'""", conn)
      
    pi_quantile = np.quantile(pi['pi'], .99)
    win_sql = f"""select win_id, chrom, start, end
                 from {poly_win_table}
                 where num_sites > 1000
                 and pi > {pi_quantile}
                 and pop = '{pop}'"""
  elif outlier_type == 'random':
    win_sql = f"""select win_id, chrom, start, end
                  from {d_win_table}
                  where num_sites > 1000
                  order by random() limit 50"""

  #populate the pandas dataframe from sqlite
  win_df = pd.read_sql(win_sql, conn)
  
  
  #get outgroup
  outgroup = pop_str.split('_')[3]


  # Intialize pop dictionary.
  idx_pop_dicc = {}
  for pop in list(set(meta_df['pop'])):
    idx_pop_dicc[pop] = meta_df['vcf_order'][meta_df['pop'] == pop]
      
  
  #define empty lists corresponding to database fields
  odwps_ids = []
  win_ids = []
  chroms = []
  positions = []
  outgroup_alleles = []
  total_sims = []
  total_sshs = []
  total_sechs = []
  total_sechs_anro = []
  total_sechs_denis = []
  total_sechs_ladigue = []
  total_sechs_marianne = []
  total_sechs_praslin = []
  total_sechs_unknown = []
  der_sims = []
  der_sshs = []
  der_sechs = []
  der_sechs_anro = []
  der_sechs_denis = []
  der_sechs_ladigue = []
  der_sechs_marianne = []
  der_sechs_praslin = []
  der_sechs_unknown = []
 
 
  
  # calculate stats for every outlier window
  odwps_id = 0
  for idx in range(win_df.shape[0]):
      win_id = win_df.win_id.values[idx]
      chrom = win_df.chrom.values[idx]
      start = win_df.start.values[idx]
      end = win_df.end.values[idx]
       
      # Extract the genotype callset and positions.
      zarr_file = zarr_prefix + '_' + chrom + '.zarr'
      callset, all_pos = load_callset_pos(chrom, zarr_file)
      # Identify the window to extract.
      wind_loc = all_pos.locate_range(start, end)
      win_gt = allel.GenotypeArray(callset[wind_loc])
      

      # Compute derived allele freq
      total_sim, der_sim = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sim'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_ssh, der_ssh = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['ssh'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_sech, der_sech = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sech'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))

      total_sech_anro, der_sech_anro = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sechanro'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_sech_denis, der_sech_denis = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sechden'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_sech_ladigue, der_sech_ladigue = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sechlad'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_sech_marianne, der_sech_marianne = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sechmari'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_sech_praslin, der_sech_praslin = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sechpras'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      total_sech_unknown, der_sech_unknown = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc['sechunk'], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))

 
      for pos in all_pos[wind_loc]:
        odwps_id += 1
        odwps_ids.append(odwps_id)
        win_ids.append(win_id)
        chroms.append(chrom)
        positions.append(pos)
        

      outgroup_alleles.extend(win_gt.take(idx_pop_dicc[outgroup], axis=1)[:, 0, 0])
      
      #add values to lists
      total_sims.extend(total_sim)
      total_sshs.extend(total_ssh)
      total_sechs.extend(total_sech)
      total_sechs_anro.extend(total_sech_anro)
      total_sechs_denis.extend(total_sech_denis)
      total_sechs_ladigue.extend(total_sech_ladigue)
      total_sechs_marianne.extend(total_sech_marianne)
      total_sechs_praslin.extend(total_sech_praslin)
      total_sechs_unknown.extend(total_sech_unknown)
      
      der_sims.extend(der_sim)
      der_sshs.extend(der_ssh)
      der_sechs.extend(der_sech)
      der_sechs_anro.extend(der_sech_anro)
      der_sechs_denis.extend(der_sech_denis)
      der_sechs_ladigue.extend(der_sech_ladigue)
      der_sechs_marianne.extend(der_sech_marianne)
      der_sechs_praslin.extend(der_sech_praslin)
      der_sechs_unknown.extend(der_sech_unknown)

  #create dataframe using lists that can be loaded into the database table           
  site_df = pd.DataFrame()    
  #add columns to site_df
  site_df["odwps_id"] = odwps_ids  
  site_df["win_id"] = win_ids  
  site_df["chrom"] = chroms  
  site_df["pos"] = positions
  site_df["mel_allele"] = outgroup_alleles
  site_df["total_alleles_sim"] = total_sims
  site_df["total_alleles_ssh"] = total_sshs
  site_df["total_alleles_sech"] = total_sechs
  site_df["total_alleles_sech_anro"] = total_sechs_anro
  site_df["total_alleles_sech_denis"] = total_sechs_denis
  site_df["total_alleles_sech_ladigue"] = total_sechs_ladigue
  site_df["total_alleles_sech_marianne"] = total_sechs_marianne
  site_df["total_alleles_sech_praslin"] = total_sechs_praslin
  site_df["total_alleles_sech_unknown"] = total_sechs_unknown
  site_df["der_alleles_sim"] = der_sims
  site_df["der_alleles_ssh"] = der_sshs
  site_df["der_alleles_sech"] = der_sechs
  site_df["der_alleles_sech_anro"] = der_sechs_anro
  site_df["der_alleles_sech_denis"] = der_sechs_denis
  site_df["der_alleles_sech_ladigue"] = der_sechs_ladigue
  site_df["der_alleles_sech_marianne"] = der_sechs_marianne
  site_df["der_alleles_sech_praslin"] = der_sechs_praslin
  site_df["der_alleles_sech_unknown"] = der_sechs_unknown

  
  
  #import site_df to db
  conn = sqlite3.connect(db_file)
  site_df.to_sql(sites_table, if_exists = 'append', index=False, con=conn)
  conn.commit()
  conn.close()
  

if __name__ == '__main__':
  main()
