import sqlite3
import os.path

#python3 verify_and_record_bams.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  mapping_dir = os.path.join(base_dir, 'mapping')
  bam_yes_file = os.path.join(mapping_dir, 'ssh_bam_yes.txt')
  bam_no_file = os.path.join(mapping_dir, 'ssh_bam_no.txt')
  
  db_file = os.path.join(base_dir, 'ssh_genotyping.db')
  bam_dir = os.path.join(mapping_dir, 'bams')
  
  conn = sqlite3.connect(db_file)    
  
  conn.execute("""drop table if exists bam_files""")
  conn.execute("""create table bam_files
                  (bf_id int primary key,
                   sample_id varchar(50),
                   bam_file varchar(255))""")
  
  conn.execute("create index idx_bf_sample_id on bam_files(sample_id)")
  
  sample_query = conn.execute(f"""select distinct sample_id
                                  from ssh_fastqs
                                  order by sample_id""")
                              
  cmd_count = 0
  bf_id = 0                              
  with open(bam_yes_file, 'w') as fy:
    with open(bam_no_file, 'w') as fn:
      for sample_id, in sample_query:
        cmd_count += 1
        
        sample_dir = os.path.join(bam_dir, sample_id)
        bam_file = os.path.join(sample_dir, sample_id + '.bam')
        bai_file = bam_file + '.bai'
        
        #check for bai files since malformed or empty bams are possible but cannot be indexed
        if os.path.isfile(bai_file):
          fy.write(f"{cmd_count}\t{bam_file}\n")
          bf_id += 1
          conn.execute(f"""insert into bam_files
                           values
                           ({bf_id}, '{sample_id}', '{bam_file}')""")                           
        else:
          fn.write(f"{cmd_count}\t{bam_file}\n")

  conn.commit()


       
if __name__ == '__main__':
  main()
