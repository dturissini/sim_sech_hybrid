import sqlite3
import os
import sys
import tempfile


#python3 load_drosophila_ref_genomes.py


def create_db_tables(db_file):
  conn = sqlite3.connect(db_file)    

  conn.execute("drop table if exists chrom_lens")
  conn.execute("""create table chrom_lens
                 (chrom varchar(50) primary key,
                  chrom_len int)""")  
  conn.close()
  
  


def get_ref_genome_chroms(db_file):
  conn = sqlite3.connect(db_file)  
    
  with tempfile.NamedTemporaryFile(mode='w') as t_drgc:
    ref_genome_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'  
    
    with open(ref_genome_file, 'r') as r:
      chrom = ''
      chrom_len = 0
      for line in r:
        if line[0] == '>':
          if chrom_len > 0:
            t_drgc.write(f"{chrom}\t{chrom_len}\n")                  
          
          chrom = line.split()[0][1:]
          chrom_len = 0
        else:
          chrom_len += len(line.strip())
          
      t_drgc.write(f"{chrom}\t{chrom_len}\n")                  
                    
    t_drgc.flush() 
    
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t_drgc.name} chrom_lens" """)
    
  conn.close()
  
  
def make_chrom_interval_lists(db_file):
  conn = sqlite3.connect(db_file)  
  ref_genome_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'  
  list_file = ref_genome_file + '.list'
  
  with open(list_file, 'w') as l:
    #test this quyery in interactive python session
    chrom_query = conn.execute(f"""select chrom
                                   from chrom_lens""")
    
    for chrom, in chrom_query:
      l.write(f"{chrom}\n")
        
  conn.close()
    

def main():
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  
  db_file = os.path.join(base_dir, 'ssh_genotyping.db')
  
  create_db_tables(db_file)
  get_ref_genome_chroms(db_file)
  make_chrom_interval_lists(db_file)
  
  



if __name__ == '__main__':
  main()
