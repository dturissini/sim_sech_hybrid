import sqlite3
import os
import glob
import re


#python3 identify_fastq_files.py




def main():
  db_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ssh_genotyping.db'
  d_win_db_file = '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
  proj_dir = '/proj/matutelb/projects/drosophila'
  htsf_dir = '/proj/matutelb/HTSF'
  fastq_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/fastqs'
  
  os.system(f"mkdir -p {fastq_dir}")
  
  conn = sqlite3.connect(db_file)
  conn_d_win = sqlite3.connect(d_win_db_file)
  
  conn.execute(f"""drop table if exists ssh_fastqs""")
  conn.execute(f"""create table ssh_fastqs
                   (sample_id varchar(50),
                    original_fastq_file varchar(255),
                    fastq_file varchar(255),
                    fastq_base varchar(255))""")
                    
  conn.execute(f"create index sshf_sample_id on ssh_fastqs(sample_id)")
  
  
  lk_species_dir = {'sim': 'simulans', 'sech': 'sechellia', 'ssh': 'sim_sech_hybrid'}
  
  
  sample_query = conn_d_win.execute(f"""select sample_id, species
                                        from sample_species
                                        where species != 'mel'
                                        and sample_id != 'SRR520350'""")
  

  for sample_id, species in sample_query:
    sample_fastq_dir = os.path.join(fastq_dir, sample_id)
    os.system(f"mkdir -p {sample_fastq_dir}")
    
    glob_str_proj = os.path.join(proj_dir, lk_species_dir[species], 'fastqs', '*' + sample_id + '_*.fastq.gz')
    glob_str_htsf = os.path.join(htsf_dir, '**', '*' + sample_id + '_*.fastq.gz')
    
    no_proj_files = True
    for original_fastq_file in glob.glob(glob_str_proj):
      no_proj_files = False
       
      fastq_file = os.path.join(sample_fastq_dir, original_fastq_file.split('/')[-1])
      fastq_base = fastq_file.split('/')[-1][:-16]
      
      os.system(f"ln -s {original_fastq_file} {fastq_file}")
      
      conn.execute(f"""insert into ssh_fastqs
                       values
                       ('{sample_id}', '{original_fastq_file}', '{fastq_file}', '{fastq_base}')""")
    
    if no_proj_files == True:
      for original_fastq_file in glob.glob(glob_str_htsf, recursive=True):
        fastq_file = os.path.join(sample_fastq_dir, original_fastq_file.split('/')[-1])
        fastq_base = fastq_file.split('/')[-1][:-16]

        os.system(f"ln -s {original_fastq_file} {fastq_file}")
        
        conn.execute(f"""insert into ssh_fastqs
                         values
                         ('{sample_id}', '{original_fastq_file}', '{fastq_file}', '{fastq_base}')""")
    

  conn.commit()

if __name__ == '__main__':
  main()
