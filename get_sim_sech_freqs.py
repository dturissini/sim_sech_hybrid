import argparse
import sqlite3
import os
import tempfile
import re
import gzip

#python3 /proj/matutelb/projects/gwas/scripts/get_sim_sech_freqs.py /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.chroms.biallelic.filtered.recode.OA.filtered.recode.vcf.gz /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_final_samplelist.tsv /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_freqs.db


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

  cur.execute("drop table if exists sim_sech_freqs")
  cur.execute("""create table sim_sech_freqs
                 (gssf_id int primary key,
                  chrom varchar(20),
                  pos int,
                  ref_allele varchar(20),
                  alt_allele varchar(20),
                  num_sim int,
                  num_sech int,
                  num_sech_denis int,
                  num_sim_sech_hybrid int,
                  sim_alt_freq decimal(5,4),
                  sech_alt_freq decimal(5,4),
                  sech_denis_alt_freq decimal(5,4),
                  sim_sech_hybrid_alt_freq decimal(5,4))""")
   
  cur.execute("create index idx_ssf_pos on sim_sech_freqs(pos)")
 
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
  gssf_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as o:
    with gzip.open(vcf_file, 'rt') as v: 
      for line in v:
        counter += 1
        if counter % 100000 == 0:
          print(counter)
        
        if line[:2] != "##":
          line = line.strip()
          values = line.split('\t')
          snp_name = values[0] + '_' + values[1]
          
          if line[0] == '#':
            fly_lines = values[9:]
          elif values[0] in lk_chrom:
            num_samples = {'simulans': 0, 'sechellia': 0, 'sech_denis': 0, 'sim_sech_hybrid': 0}
            num_alts = {'simulans': 0, 'sechellia': 0, 'sech_denis': 0, 'sim_sech_hybrid': 0}
            for i, sample_str in enumerate(values[9:]):     
              alt_count = 0
              attributes = sample_str.split(':')
              if attributes[2] == '.':
                attributes[2] = 0
                
              if int(attributes[2]) > min_depth:
                for allele in re.split(r"/|\|", attributes[0]):
                  if allele == '1':
                    alt_count += 1
                
                num_samples[lk_species[fly_lines[i]]] += 1
                num_alts[lk_species[fly_lines[i]]] += alt_count
                if lk_species[fly_lines[i]] == 'sechellia' and 'Denis' in fly_lines[i]:
                  num_samples['sech_denis'] += 1
                  num_alts['sech_denis'] += alt_count
                        
            sim_alt_freq = -1
            sech_alt_freq = -1
            sech_denis_alt_freq = -1
            sim_sech_hybrid_alt_freq = -1
            if num_samples['simulans'] > 0: 
              sim_alt_freq = num_alts['simulans'] / num_samples['simulans'] / 2
                        
            if num_samples['sechellia'] > 0: 
              sech_alt_freq = num_alts['sechellia'] / num_samples['sechellia'] / 2

            if num_samples['sech_denis'] > 0: 
              sech_denis_alt_freq = num_alts['sech_denis'] / num_samples['sech_denis'] / 2

            if num_samples['sim_sech_hybrid'] > 0: 
              sim_sech_hybrid_alt_freq = num_alts['sim_sech_hybrid'] / num_samples['sim_sech_hybrid'] / 2
                            
              gssf_id += 1        
              o.write(f"{gssf_id}\t{lk_chrom[values[0]]}\t{values[1]}\t{values[3]}\t{values[4]}\t{num_samples['simulans']}\t{num_samples['sechellia']}\t{num_samples['sech_denis']}\t{num_samples['sim_sech_hybrid']}\t{sim_alt_freq}\t{sech_alt_freq}\t{sech_denis_alt_freq}\t{sim_sech_hybrid_alt_freq}\n")                  
    
    o.flush() 
    
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {o.name} sim_sech_freqs" """)


  conn = sqlite3.connect(db_file)    
  cur = conn.cursor()                            

  cur.execute("""update sim_sech_freqs 
                 set sim_alt_freq = null
                 where sim_alt_freq = -1""")

  cur.execute("""update sim_sech_freqs 
                 set sech_alt_freq = null
                 where sech_alt_freq = -1""")

  cur.execute("""update sim_sech_freqs 
                 set sech_denis_alt_freq = null
                 where sech_denis_alt_freq = -1""")

  cur.execute("""update sim_sech_freqs 
                 set sim_sech_hybrid_alt_freq = null
                 where sim_sech_hybrid_alt_freq = -1""")
  conn.commit()
  conn.close()  

if __name__ == '__main__':
  main()
