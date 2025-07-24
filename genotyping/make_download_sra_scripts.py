import sqlite3
import os
import sys


#python3 make_download_sra_scripts.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  sra_dir = os.path.join(base_dir, 'sra')
  sample_sh_dir = os.path.join(sra_dir, 'sample_sra_scripts')
  fastq_dir = os.path.join(base_dir, 'fastqs')
  log_dir = os.path.join(sra_dir, 'logs')
  sra_tmp_dir = os.path.join(sra_dir, 'sra_tmp')

  db_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ssh_genotyping.db'
  master_sh_file = os.path.join(sra_dir, 'ssh_srrs_fasterq.sh')
  sbatch_file = os.path.join(sra_dir, 'ssh_srrs_fasterq.sbatch')
 
  os.system(f"mkdir -p {sample_sh_dir}")
  os.system(f"mkdir -p {fastq_dir}")
  os.system(f"mkdir -p {log_dir}")
  os.system(f"mkdir -p {sra_tmp_dir}")
  
  conn = sqlite3.connect(db_file)   
  conn.execute(f"""attach database '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db' as d""") 
  
  missing_fastq_query = conn.execute(f"""select species, sample_id
                                         from sample_species
                                         where sample_id not in (select distinct sample_id from ssh_fastqs)
                                         and substr(sample_id, 1, 3) in ('ERR', 'SRR')""")

  total_cmds = 0
  with open(master_sh_file, 'w') as m:
    for species, sample_id in missing_fastq_query:
      total_cmds += 1
      
      sample_dir = os.path.join(fastq_dir, sample_id)
      sample_sh_file = os.path.join(sample_sh_dir, sample_id + '_get_fastq.sh')
      outfile_base = os.path.join(sample_dir, sample_id)

      os.system(f"mkdir -p {sample_dir}")
      
      m.write(f"{sample_sh_file}\n")
      
      with open(sample_sh_file, 'w') as i:       
        outfile = outfile_base + '.fastq'

        i.write(f"fasterq-dump -o {outfile} -e 4 --split-files -f {sample_id}\n")        
        i.write(f"\ngzip {outfile_base}*.fastq\n")
        
        conn.execute(f"""insert into ssh_fastqs
                         values
                         ('{sample_id}', '{outfile_base}_1.fastq.gz', '{outfile_base}_1.fastq.gz', '{sample_id}')""")
        
        conn.execute(f"""insert into ssh_fastqs
                         values
                         ('{sample_id}', '{outfile_base}_2.fastq.gz', '{outfile_base}_2.fastq.gz', '{sample_id}')""")
        

      os.system(f"chmod +x {sample_sh_file}")

  conn.commit()



  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J ssh_sra\n")
    o.write(f"#SBATCH --cpus-per-task=4\n")
    o.write(f"#SBATCH -t 8:00:00\n")
    o.write(f"#SBATCH --mem=4g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}%50\n")
    o.write(f"#SBATCH -o {log_dir}/ssh_srrs_fasterq_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/ssh_srrs_fasterq_%A_%a.err\n\n")
    
    o.write(f"ml sratoolkit/3.1.1\n\n")

    o.write(f"fasterq_cmd=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {master_sh_file})\n")
    o.write(f"""eval "${{fasterq_cmd}}"\n""")
  
       
if __name__ == '__main__':
  main()
