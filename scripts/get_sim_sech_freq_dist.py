import argparse
import sqlite3
import os
import tempfile
import re
import gzip
import math
from datetime import datetime

#python3 /proj/matutelb/projects/gwas/scripts/get_sim_sech_freq_dist.py /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.vcf.gz /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_final_samplelist.tsv /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_freq_dist.db


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

  cutoff_sim = 60
  cutoff_sech = 20
  cutoff_ssh = 10

  cur.execute("drop table if exists sim_sech_freq_dist")
  cur.execute("""create table sim_sech_freq_dist
                 (ssfd_id int primary key,
                  sim_alt_freq_bin decimal(5,4),
                  sech_alt_freq_bin decimal(5,4),
                  sim_sech_hybrid_alt_freq_bin decimal(5,4),
                  num_sites int)""")
    
  conn.close()

  min_depth = 6
  lk_chrom = {'NC_052520.2': '2L',
              'NC_052521.2': '2R',
              'NC_052522.2': '3L',
              'NC_052523.2': '3R',
              'NC_052524.2': '4',
              'NC_052525.2': 'X'}
 
 
#process species_list_file
  lk_species = {}
  with open(species_list_file, 'r') as l: 
    for line in l:
      line = line.strip()
      values = line.split('\t')
      lk_species[values[0]] = values[1]


#initiate freq_bins   
  freq_bins = {}
    
  num_bins = 20
  bins = list(x * 1 / num_bins for x in range(0, num_bins + 1))
  for sim_bin in bins:
    freq_bins[sim_bin] = {}
    for sech_bin in bins:
      freq_bins[sim_bin][sech_bin] = {}
      for ssh_bin in bins:
        freq_bins[sim_bin][sech_bin][ssh_bin] = 0
  

  
  counter = 0
  with tempfile.NamedTemporaryFile(mode='w') as o:
    with gzip.open(vcf_file, 'rt') as v: 
      for line in v:
        counter += 1
        if counter % 100000 == 0:
          print(counter, datetime.now())
        
        if line[:2] != "##":
          line = line.strip()
          values = line.split('\t')
          snp_name = values[0] + '_' + values[1]
          
          if line[0] == '#':
            fly_lines = values[9:]
          elif values[0] in lk_chrom:
            num_samples = {'simulans': 0, 'sechellia': 0, 'sim_sech_hybrid': 0}
            num_alts = {'simulans': 0, 'sechellia': 0, 'sim_sech_hybrid': 0}
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
                        
            if num_samples['simulans'] > cutoff_sim and num_samples['sechellia'] > cutoff_sech and num_samples['sim_sech_hybrid'] > cutoff_ssh:
              sim_alt_freq = num_alts['simulans'] / num_samples['simulans'] / 2
              sech_alt_freq = num_alts['sechellia'] / num_samples['sechellia'] / 2
              sim_sech_hybrid_alt_freq = num_alts['sim_sech_hybrid'] / num_samples['sim_sech_hybrid'] / 2
              
              sim_bin = get_bin(sim_alt_freq, num_bins)
              sech_bin = get_bin(sech_alt_freq, num_bins)
              ssh_bin = get_bin(sim_sech_hybrid_alt_freq, num_bins)
              
              freq_bins[sim_bin][sech_bin][ssh_bin] += 1
              
              
                            
    ssfd_id = 0
    for sim_bin in freq_bins:
      for sech_bin in freq_bins[sim_bin]:
        for ssh_bin in freq_bins[sim_bin][sech_bin]:    
          ssfd_id += 1
          o.write(f"{ssfd_id}\t{sim_bin}\t{sech_bin}\t{ssh_bin}\t{freq_bins[sim_bin][sech_bin][ssh_bin]}\n")                  
    
    o.flush() 
    
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {o.name} sim_sech_freq_dist" """)



def get_bin(freq, num_bins):
  bin = math.floor(freq * num_bins) / num_bins    
  return bin

if __name__ == '__main__':
  main()
