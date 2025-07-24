import sqlite3
import os
import sys


#python3 make_gvcf_scripts.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'

  mapping_dir = os.path.join(base_dir, 'mapping')
  bam_dir = os.path.join(mapping_dir, 'bams')
  gvcf_dir = os.path.join(base_dir, 'gvcfs')
  gvcf_sh_dir = os.path.join(base_dir, 'gvcf_sh')
  log_dir = os.path.join(base_dir, 'logs')

  ref_genome_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'  
  master_sh_file = os.path.join(base_dir, 'make_gvcfs.sh')
  sbatch_file = os.path.join(base_dir, 'make_gvcfs.sbatch')
 
  os.system(f"mkdir -p {gvcf_dir}")
  os.system(f"mkdir -p {gvcf_sh_dir}")
  os.system(f"mkdir -p {log_dir}")
  
  db_file = os.path.join(base_dir, 'ssh_genotyping.db')
  conn = sqlite3.connect(db_file)    
  sample_query = conn.execute(f"""select distinct sample_id, bam_file
                                  from bam_files
                                  order by sample_id""")
  
                              
  total_cmds = 0
  with open(master_sh_file, 'w') as m:
    for sample_id, bam_file in sample_query:
      total_cmds += 1
      
      gvcf_sample_dir = os.path.join(gvcf_dir, sample_id)
      gvcf_file = os.path.join(gvcf_sample_dir, sample_id + '.gvcf.gz')

      ind_sh_file = os.path.join(gvcf_sh_dir, sample_id + '_make_gvcf.sh')

      os.system(f"mkdir -p {gvcf_sample_dir}")
      
      m.write(f"{ind_sh_file}\n")
      
      with open(ind_sh_file, 'w') as i:
        
        i.write(f"set -o pipefail\n\n")
        i.write(f"""gatk --java-options "-Xmx8g" HaplotypeCaller  \\\n""")
        i.write(f"""-R {ref_genome_file} \\\n""")
        i.write(f"""-I {bam_file} \\\n""")
        i.write(f"""-O {gvcf_file} \\\n""")
        i.write(f"""-ERC GVCF\n""")
          
      os.system(f"chmod +x {ind_sh_file}")



  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J gvcf\n")
    o.write(f"#SBATCH -t 240:00:00\n")
    o.write(f"#SBATCH --mem=8g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}%100\n")
    o.write(f"#SBATCH -o {log_dir}/gvcf_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/gvcf_%A_%a.err\n\n")
    
    o.write(f"ml gatk/4.6.1.0\n")

    o.write(f"gatk_sh=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {master_sh_file})\n")
    o.write(f"""eval "${{gatk_sh}}"\n""")


       
if __name__ == '__main__':
  main()
