import os
import sys
import glob
import tempfile
import sqlite3
import pandas as pd
import re
from datetime import datetime


#python3 get_derived_alleles_from_read_mafs.py 



def main():
  per_der_db = 'inversion_check_per_der.db'  

  conn = sqlite3.connect(per_der_db)
  
  cur = conn.cursor()                            

  cur.execute("drop table if exists inversion_check_per_der_reads")
  cur.execute("""create table inversion_check_per_der_reads
                 (icpdt_id int primary key,
                  sample_id varcahr(50),
                  inv_name varchar(20),
                  read_name varchar(20),
                  num_der_sites int,
                  num_anc_sites int,
                  total_sites int,
                  num_anro_b3_alleles int)""")
   
  cur.execute("create index idx_icpdt_sample_id on inversion_check_per_der_reads(sample_id)")
  cur.execute("create index idx_icpdt_inv_name on inversion_check_per_der_reads(inv_name)")


  cur.execute("drop table if exists inversion_check_read_alleles")
  cur.execute("""create table inversion_check_read_alleles
                 (icra_id int primary key,
                  sample_id varcahr(50),
                  inv_name varchar(20),
                  read_name varchar(20),
                  chrom varchar(5),
                  pos int,
                  allele varchar(1),
                  der_site int,
                  anc_site int,
                  anro_site int)""")
   
  cur.execute("create index idx_icra_sample_id on inversion_check_read_alleles(sample_id)")
  cur.execute("create index idx_icra_inv_name on inversion_check_read_alleles(inv_name)")
  cur.execute("create index idx_icra_pos on inversion_check_read_alleles(pos)")

    
  sites_alleles_df = pd.read_sql(f"""select distinct chrom, pos, outgroup_allele, der_allele, anro_allele
                                     from inversion_sech_alleles""", conn)
  
  conn.close()
  
  maf_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly/inversion_check/mafs_reads'
  
  icpdt_id = 0
  icra_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t_sites:
    with tempfile.NamedTemporaryFile(mode='w') as t_reads:
      for maf_file in glob.glob(os.path.join(maf_dir, '*.maf')):
        if icpdt_id % 1000 == 0:
          print(icpdt_id, datetime.now())
          
        num_der_sites = 0
        num_anc_sites = 0
        num_anro_b3_alleles = 0
        with open(maf_file, 'r') as m:
          m_maf = re.search(r'.+/mafs_reads/(.+_\d+M)_(.+)_([a-zA-Z0-9\-]+).maf', maf_file)
          inv_name, sample_id, read_name = m_maf.groups()

          for line in m:
            if line[0] == 'a':
              ref_line = next(m)
              read_line = next(m)
              
              ref_s, ref_chrom, ref_start, ref_aln_len, ref_strand, ref_chrom_len, ref_seq = ref_line.strip().split()
              read_s, read, read_start, read_aln_len, read_strand, read_len, read_seq = read_line.strip().split()
              
              ref_start = int(ref_start) + 1
              read_start = int(read_start) + 1
              ref_aln_len = int(ref_aln_len)
              
              aln_sites_df = sites_alleles_df[(sites_alleles_df['chrom'] == ref_chrom) & (sites_alleles_df['pos'] >= ref_start) & (sites_alleles_df['pos'] <= ref_start + ref_aln_len)] 
              
              for index, row in aln_sites_df.iterrows():
                aln_pos = row['pos'] - ref_start
                ref_allele = ref_seq[aln_pos:(aln_pos + 1)]
                read_allele = read_seq[aln_pos:(aln_pos + 1)]
                
                der_site = 0
                anc_site = 0
                anro_site = 0                
                if read_allele == row['der_allele']:
                  num_der_sites += 1
                  der_site = 1
                elif read_allele == row['outgroup_allele']:
                  num_anc_sites += 1
                  anc_site = 1
    
                if read_allele == row['anro_allele']:
                  num_anro_b3_alleles += 1
                  anro_site = 1

                icra_id += 1
                t_sites.write(f"{icra_id}\t{sample_id}\t{inv_name}\t{read_name}\t{row['chrom']}\t{row['pos']}\t{read_allele}\t{der_site}\t{anc_site}\t{anro_site}\n")                

                
          icpdt_id += 1
          total_sites = num_der_sites + num_anc_sites
          
          t_reads.write(f"{icpdt_id}\t{sample_id}\t{inv_name}\t{read_name}\t{num_der_sites}\t{num_anc_sites}\t{total_sites}\t{num_anro_b3_alleles}\n")                
                             
      t_reads.flush()      
      os.system(f"""sqlite3 {per_der_db} ".mode tabs" ".import {t_reads.name} inversion_check_per_der_reads" """)

    t_sites.flush()      
    os.system(f"""sqlite3 {per_der_db} ".mode tabs" ".import {t_sites.name} inversion_check_read_alleles" """)

if __name__ == '__main__':
  main()
