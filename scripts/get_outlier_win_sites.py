# Import packages.
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
import sqlite3


#create d stats table
def create_windows_table(conn, win_size, outlier_type, sites_table, sites_alleles_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {sites_table}")
  cur.execute(f"""create table {sites_table}
                  (odwps_id int primary key,
                   win_id varchar(30),
                   chrom varchar(20),
                   pos int,
                   outgroup_allele_numeric varchar(1),
                   pop varchar(50),
                   total_alleles int,
                   der_alleles int)""")

                   
  cur.execute(f"create index idx_odwps_pos_{win_size}{outlier_type} on {sites_table}(pos)")
  cur.execute(f"create index idx_odwps_win_id_{win_size}{outlier_type} on {sites_table}(win_id)")
  cur.execute(f"create index idx_odwps_pop_{win_size}{outlier_type} on {sites_table}(pop)")
  
  
  cur.execute(f"drop table if exists {sites_alleles_table}")
  cur.execute(f"""create table {sites_alleles_table}
                  (odwpsa_id int primary key,
                   chrom varchar(20),
                   pos int,
                   outgroup_allele varchar(1),
                   der_allele varchar(1))""")

                   
  cur.execute(f"create index idx_odwpsa_pos_{win_size}{outlier_type} on {sites_alleles_table}(pos)")
  
  cur.close()
  


# Define a function to load genotyope and positions arrays.
def load_callset_pos(chrom, zarr_file):
    # Load the vcf file.
    callset = zarr.open_group(zarr_file, mode='r')
    # Extract the genotypes.
    geno = callset[f'{chrom}/calldata/GT']
    # Extract the ref and alt alleles.
    refs = callset[f'{chrom}/variants/REF']
    alts = callset[f'{chrom}/variants/ALT']
    # Load the positions.
    pos = allel.SortedIndex(callset[f'{chrom}/variants/POS'])
    return geno, pos, refs, alts


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
  sites_alleles_table = "outlier_" + outlier_type + "_win_sites_alleles_" + str(win_size) + '_' + pop_str

  create_windows_table(conn, win_size, outlier_type, sites_table, sites_alleles_table)  


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
    pi_pop = outlier_type.split('_')[1]
    pi = pd.read_sql(f"""select pi
                         from {poly_win_table}
                         where num_sites > 1000
                         and pop = '{pi_pop}'""", conn)
      
    pi_quantile = np.quantile(pi['pi'], .99)
    win_sql = f"""select win_id, chrom, start, end
                 from {poly_win_table}
                 where num_sites > 1000
                 and pi > {pi_quantile}
                 and pop = '{pi_pop}'"""
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
  for pop_i in list(set(meta_df['pop'])):
    idx_pop_dicc[pop_i] = meta_df['vcf_order'][meta_df['pop'] == pop_i]
      
  focal_pops = list(set(meta_df['pop']))
  focal_pops.remove(outgroup)
  
  

  # calculate stats for every outlier window
  odwps_id = 0
  odwpsa_id = 0
  for idx in range(win_df.shape[0]):
    win_id = win_df.win_id.values[idx]
    chrom = win_df.chrom.values[idx]
    start = win_df.start.values[idx]
    end = win_df.end.values[idx]
     
    # Extract the genotype callset and positions.
    zarr_file = zarr_prefix + '_' + chrom + '.zarr'
    callset, all_pos, refs, alts = load_callset_pos(chrom, zarr_file)
    # Identify the window to extract.
    wind_loc = all_pos.locate_range(start, end)
    win_gt = allel.GenotypeArray(callset[wind_loc])
    
    for pop_i in focal_pops:
      #define empty lists corresponding to database fields
      odwps_ids = []
      win_ids = []
      chroms = []
      positions = []
      outgroup_alleles_numeric = []
      pops = []
      total_alleles = []
      der_alleles = []

    
      # Compute derived allele freq
      total_alleles, der_alleles = get_der_allele_counts(gt=win_gt.take(idx_pop_dicc[pop_i], axis=1), outgroup_gt=win_gt.take(idx_pop_dicc[outgroup], axis=1))
      
      #get outgroup alleles
      outgroup_alleles_numeric.extend(win_gt.take(idx_pop_dicc[outgroup], axis=1)[:, 0, 0])
      
      for pos in all_pos[wind_loc]:
        odwps_id += 1
        odwps_ids.append(odwps_id)
        win_ids.append(win_id)
        chroms.append(chrom)
        positions.append(pos)
        pops.append(pop_i)
        
      #create dataframe using lists that can be loaded into the database table           
      site_df = pd.DataFrame()    

      #add columns to site_df      
      site_df["odwps_id"] = odwps_ids  
      site_df["win_id"] = win_ids  
      site_df["chrom"] = chroms  
      site_df["pos"] = positions
      site_df["outgroup_allele_numeric"] = outgroup_alleles_numeric
      site_df["pop"] = pops
      site_df["total_alleles"] = total_alleles
      site_df["der_alleles"] = der_alleles
      
      
           
      #import site_df db
      conn = sqlite3.connect(db_file)
      site_df.to_sql(sites_table, if_exists = 'append', index=False, con=conn)

      if pop_i == focal_pops[0]:
        #get outgroup and derived alleles
        odwpsa_ids = []
        allele_chroms = []
        allele_positions = []
        outgroup_alleles = []
        der_alleles = []
        
        outgroup_allele_i = 0
        for i in range(wind_loc.start, wind_loc.stop):
          alleles = [refs[i][0], alts[i][0]]        
          der_allele_numeric = abs(outgroup_alleles_numeric[outgroup_allele_i] - 1)
          outgroup_allele = alleles[outgroup_alleles_numeric[outgroup_allele_i]]
          der_allele = alleles[der_allele_numeric]
          
          outgroup_allele_i += 1
          odwpsa_id += 1
        
          odwpsa_ids.append(odwpsa_id)
          allele_chroms.append(chrom)
          allele_positions.append(all_pos[i])
          
          outgroup_alleles.append(outgroup_allele)
          der_alleles.append(der_allele)
          
        
             
        #create dataframe using lists that can be loaded into the database table           
        site_allele_df = pd.DataFrame()    
        
        #add columns to site_df      
        site_allele_df["odwpsa_id"] = odwpsa_ids  
        site_allele_df["chrom"] = allele_chroms  
        site_allele_df["pos"] = allele_positions
        site_allele_df["outgroup_allele"] = outgroup_alleles
        site_allele_df["der_allele"] = der_alleles
        
        #import site_df_allele_to db
        site_allele_df.to_sql(sites_alleles_table, if_exists = 'append', index=False, con=conn)
    
  conn.commit()
  conn.close()
  

if __name__ == '__main__':
  main()
