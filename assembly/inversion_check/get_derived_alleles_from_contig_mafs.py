import os
import sys
import glob
import tempfile
import sqlite3
import pandas as pd
import re


#python3 get_derived_alleles_from_contig_mafs.py 



def main():
  per_der_db = 'inversion_check_per_der.db'  

  conn = sqlite3.connect(per_der_db)  
  
  cur = conn.cursor()                            

  cur.execute("drop table if exists inversion_check_per_der")
  cur.execute("""create table inversion_check_per_der
                 (icpd_id int primary key,
                  sample_id varcahr(50),
                  inv_name varchar(20),
                  contig_name varchar(20),
                  num_der_sites int,
                  num_anc_sites int,
                  total_sites int,
                  num_anro_b3_alleles int)""")
   
  cur.execute("create index idx_icpd_sample_id on inversion_check_per_der(sample_id)")
  cur.execute("create index idx_icpd_inv_name on inversion_check_per_der(inv_name)")
  
  sites_alleles_df = pd.read_sql(f"""select distinct chrom, pos, outgroup_allele, der_allele, anro_allele
                                     from inversion_sech_alleles""", conn)
  
  conn.close()
  
  maf_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly/inversion_check/mafs'
  
  icpd_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t:
    for maf_file in glob.glob(os.path.join(maf_dir, '*.maf')):
      print(maf_file)
      num_der_sites = 0
      num_anc_sites = 0
      num_anro_b3_alleles = 0
      with open(maf_file, 'r') as m:
        m_maf = re.search(r'.+/mafs/(.._\d+M)_(.+)_([a-zA-Z0-9]+).maf', maf_file)
        inv_name, sample_id, contig_name = m_maf.groups()
        
        for line in m:
          if line[0] == 'a':
            ref_line = next(m)
            contig_line = next(m)
            
            ref_s, ref_chrom, ref_start, ref_aln_len, ref_strand, ref_chrom_len, ref_seq = ref_line.strip().split()
            contig_s, contig, contig_start, contig_aln_len, contig_strand, contig_len, contig_seq = contig_line.strip().split()
            
            ref_start = int(ref_start) + 1
            contig_start = int(contig_start) + 1
            ref_aln_len = int(ref_aln_len)
            
            aln_sites_df = sites_alleles_df[(sites_alleles_df['chrom'] == ref_chrom) & (sites_alleles_df['pos'] >= ref_start) & (sites_alleles_df['pos'] <= ref_start + ref_aln_len)] 
            
            for index, row in aln_sites_df.iterrows():
              aln_pos = row['pos'] - ref_start
              ref_allele = ref_seq[aln_pos:(aln_pos + 1)]
              contig_allele = contig_seq[aln_pos:(aln_pos + 1)]
              
              
              if contig_allele == row['der_allele']:
                num_der_sites += 1
              elif contig_allele == row['outgroup_allele']:
                num_anc_sites += 1

              if contig_allele == row['anro_allele']:
                num_anro_b3_alleles += 1

        
        icpd_id += 1
        total_sites = num_der_sites + num_anc_sites
        
        t.write(f"{icpd_id}\t{sample_id}\t{inv_name}\t{contig_name}\t{num_der_sites}\t{num_anc_sites}\t{total_sites}\t{num_anro_b3_alleles}\n")                
                           
    t.flush()      
    os.system(f"""sqlite3 {per_der_db} ".mode tabs" ".import {t.name} inversion_check_per_der" """)

if __name__ == '__main__':
  main()
