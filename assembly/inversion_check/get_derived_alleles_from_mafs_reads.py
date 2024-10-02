import os
import sys
import glob
import tempfile
import sqlite3
import pandas as pd
import re


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
  
##  sites_alleles_df = pd.read_sql(f"""select sa.chrom, sa.pos, outgroup_allele, der_allele, a.num_der_alleles num_anro_b3_der_alleles
##                                     from outlier_pi_sech_win_sites_alleles_50000_sim_ssh_sech_mel sa, outlier_pi_sech_win_alleles_sech_50000_sim_ssh_sech_mel a
##                                     where sa.chrom = a.chrom
##                                     and sa.pos = a.pos
##                                     and a.num_der_alleles != -2
##                                     and a.sample_id = 'SECH_3-sech_Anro_B3_TTAGGC_L001'
##                                     and exists (select 'x' from outlier_pi_sech_win_sites_50000_sim_ssh_sech_mel s
##                                                 where sa.chrom = s.chrom
##                                                 and sa.pos = s.pos
##                                                 and 1.0 * der_alleles / total_alleles > 0
##                                                 and 1.0 * der_alleles / total_alleles < 1 
##                                                 and pop = 'sech')""", conn)
  
  sites_alleles_df = pd.read_sql(f"""select distinct chrom, pos, outgroup_allele, der_allele, anro_allele
                                     from inversion_sech_alleles""", conn)
  
  conn.close()
  
  maf_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly/inversion_check/mafs_reads'
  
  icpdt_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t:
    for maf_file in glob.glob(os.path.join(maf_dir, '*.maf')):
      print(maf_file)
      num_der_sites = 0
      num_anc_sites = 0
      num_anro_b3_alleles = 0
      with open(maf_file, 'r') as m:
        m_maf = re.search(r'.+/mafs_reads/(.._\d+M)_(.+)_([a-zA-Z0-9\-]+).maf', maf_file)
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
              
              
##              if read_allele == row['der_allele']:
##                num_der_sites += 1
##                if row['num_anro_b3_der_alleles'] == 2:
##                  num_anro_b3_alleles += 1
##              elif read_allele == row['outgroup_allele']:
##                num_anc_sites += 1
##                if row['num_anro_b3_der_alleles'] == 0:
##                  num_anro_b3_alleles += 1

              if read_allele == row['der_allele']:
                num_der_sites += 1
              elif read_allele == row['outgroup_allele']:
                num_anc_sites += 1

              if read_allele == row['anro_allele']:
                num_anro_b3_alleles += 1
        
        icpdt_id += 1
        total_sites = num_der_sites + num_anc_sites
        
        t.write(f"{icpdt_id}\t{sample_id}\t{inv_name}\t{read_name}\t{num_der_sites}\t{num_anc_sites}\t{total_sites}\t{num_anro_b3_alleles}\n")                
                           
    t.flush()      
    os.system(f"""sqlite3 {per_der_db} ".mode tabs" ".import {t.name} inversion_check_per_der_reads" """)

if __name__ == '__main__':
  main()
