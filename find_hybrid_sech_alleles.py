import argparse
import sqlite3
import os
import tempfile
import re
import gzip
from datetime import datetime

#python3 /proj/matutelb/projects/gwas/scripts/find_hybrid_sech_alleles.py /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.chroms.biallelic.filtered.recode.OA.filtered.recode.vcf.gz /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_final_samplelist.tsv /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_diff_sites.db


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_file")
    parser.add_argument("species_list_file")
    parser.add_argument("db_file")

    return parser.parse_args()


def main():
  args = parse_args()
  vcf_file = args.vcf_file  
  species_list_file = args.species_list_file  
  db_file = args.db_file  

  conn = sqlite3.connect(db_file)
    
  cur = conn.cursor()                            

  cur.execute("drop table if exists hybrid_sech_alleles")
  cur.execute("""create table hybrid_sech_alleles
                 (hsa_id int primary key,
                  chrom varchar(20),
                  pos int,
                  hybrid_line varchar(60),
                  num_sech_alleles int)""")
   
  cur.execute("create index idx_hsa_pos on hybrid_sech_alleles(pos)")
 
  diff_sites = set()
  diff_site_query = cur.execute("""select chrom, pos 
                                   from sim_sech_diff_sites
                                   where sim_ref_freq = 1
                                   and sech_ref_freq = 0
                                   and num_sech >= 10
                                   and num_sim >= 40""")
                                   
  for (chrom, pos) in diff_site_query:
    snp = chrom + '_' + str(pos)
    diff_sites.add(snp)                                 
 
  conn.close()

  min_depth = 6
  lk_chrom = {'1': '2L',
              '2': '2R',
              '3': '3L',
              '4': '3R',
              '5': '4',
              '6': 'X'}
 


#process species_list_file
  lk_species = {}
  with open(species_list_file, 'r') as l: 
    for line in l:
      line = line.strip()
      values = line.split('\t')
      lk_species[values[0]] = values[1]

  
  counter = 0
  hsa_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as o:
    with gzip.open(vcf_file, 'rt') as v: 
      for line in v:
        if counter % 100000 == 0:
          print(counter, datetime.now())
        counter += 1
        
        if line[:1] != "##":
          line = line.strip()
          values = line.split('\t')
          
          if line[0] == '#':
            fly_lines = values[9:]
          elif values[0] in lk_chrom:
            snp = lk_chrom[values[0]] + '_' + values[1]
            if snp in diff_sites:              
              for i, sample_str in enumerate(values[9:]):     
                attributes = sample_str.split(':')
                if attributes[2] == '.':
                  attributes[2] = 0
                
                if int(attributes[2]) > min_depth and lk_species[fly_lines[i]] == 'sim_sech_hybrid':
                  num_sech_alleles = 0
                  for allele in re.split(r"/|\|", attributes[0]):
                    if allele == '1':
                      num_sech_alleles += 1
                  
                  if num_sech_alleles > 0:
                    hsa_id += 1        
                    o.write(f"{hsa_id}\t{lk_chrom[values[0]]}\t{values[1]}\t{fly_lines[i]}\t{num_sech_alleles}\n")                  
                        
    o.flush() 
    
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {o.name} hybrid_sech_alleles" """)

   

  

if __name__ == '__main__':
  main()
