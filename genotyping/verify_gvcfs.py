import sqlite3
import os.path

#python3 verify_raw_gvcfs.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  gvcf_dir = os.path.join(base_dir, 'gvcfs')
  gvcf_yes_file = os.path.join(gvcf_dir, 'ssh_gvcf_yes.txt')
  gvcf_no_file = os.path.join(gvcf_dir, 'ssh_gvcf_no.txt')
  
  db_file = os.path.join(base_dir, 'ssh_genotyping.db')

  
  conn = sqlite3.connect(db_file)    
    
  sample_query = conn.execute(f"""select distinct sample_id
                                  from bam_files
                                  order by sample_id""")
                              
  cmd_count = 0
  with open(gvcf_yes_file, 'w') as fy:
    with open(gvcf_no_file, 'w') as fn:
      for sample_id, in sample_query:
        cmd_count += 1
        
        sample_dir = os.path.join(gvcf_dir, sample_id)
        gvcf_file = os.path.join(sample_dir, sample_id + '.gvcf.gz')
        tbi_file = gvcf_file + '.tbi'
        
        #check for tbi files since malformed gvcfs cannot be indexed
        if os.path.isfile(tbi_file):
          fy.write(f"{cmd_count}\t{gvcf_file}\n")
        else:
          fn.write(f"{cmd_count}\t{gvcf_file}\n")

  conn.commit()

       
if __name__ == '__main__':
  main()
