import sqlite3
import os
import sys


#python3 make_symlink_and_index_bam_scripts.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  original_bam_dir = '/proj/matutelb/users/stuckert/simulans/bams'
  
  
  bam_dir = os.path.join(base_dir, 'bams')
  bam_sh_dir = os.path.join(base_dir, 'bam_sh')
  log_dir = os.path.join(base_dir, 'logs')

  master_sh_file = os.path.join(base_dir, 'symlink_and_index_bams.sh')
  sbatch_file = os.path.join(base_dir, 'symlink_and_index_bams.sbatch')
 
  os.system(f"mkdir -p {bam_dir}")
  os.system(f"mkdir -p {bam_sh_dir}")
  os.system(f"mkdir -p {log_dir}")
  
  db_file = '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
  conn = sqlite3.connect(db_file)    
  sample_query = conn.execute(f"""select sample_id
                                  from sample_species
                                  where species != 'mel'""")
  
                              
  total_cmds = 0
  with open(master_sh_file, 'w') as m:
    for sample_id, in sample_query:
      total_cmds += 1
      
      original_bam_file = os.path.join(original_bam_dir, sample_id + '.dedupd.bam')    

      bam_sample_dir = os.path.join(bam_dir, sample_id)
      bam_file = os.path.join(bam_sample_dir, sample_id + '.dedupd.bam')    
        
      ind_sh_file = os.path.join(bam_sh_dir, sample_id + '_symlink_and_index_bam.sh')

      os.system(f"mkdir -p {bam_sample_dir}")
      
      m.write(f"{ind_sh_file}\n")
      
      with open(ind_sh_file, 'w') as i:        
        i.write(f"set -o pipefail\n\n")
        i.write(f"ln -s {original_bam_file} {bam_file}\n")
        i.write(f"samtools index {bam_file}\n")
                  
      os.system(f"chmod +x {ind_sh_file}")



  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J symlink_index_bams\n")
    o.write(f"#SBATCH -t 4:00:00\n")
    o.write(f"#SBATCH --mem=4g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}%100\n")
    o.write(f"#SBATCH -o {log_dir}/symlink_index_bam_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/symlink_index_bam_%A_%a.err\n\n")
    
    o.write(f"ml samtools/1.21\n\n")

    o.write(f"bam_sh=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {master_sh_file})\n")
    o.write(f"""eval "${{bam_sh}}"\n""")


       
if __name__ == '__main__':
  main()
