import os
import sys


#python3 make_filtered_vcf_by_chrom_script.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  vcf_dir = os.path.join(base_dir, 'vcfs')
  log_dir = os.path.join(base_dir, 'logs')  
  filtered_vcf_dir = os.path.join(base_dir, 'filtered_vcfs')
  
  sbatch_file = os.path.join(base_dir, 'add_mel_to_chrom_vcfs.sbatch')

  
  lk_chrom = ['2L',
              '2R',
              '3L',
              '3R',
              'X',
              '4']
           
  total_cmds = len(lk_chrom)
  
  chrom_str = ''
  for chrom in lk_chrom:     
    chrom_str += '"' + chrom + '" '

  chrom_str = chrom_str.rstrip()
  
                              
  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J add_mel_to_chrom_vcfs\n")
    o.write(f"#SBATCH -t 96:00:00\n")
    o.write(f"#SBATCH --mem=8g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}\n")
    o.write(f"#SBATCH -o {log_dir}/add_mel_to_chrom_vcfs_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/add_mel_to_chrom_vcfs_%A_%a.err\n\n")
    
    
    o.write(f"""chroms=({chrom_str})\n\n""")

    o.write(f"i=$((${{SLURM_ARRAY_TASK_ID}} - 1))\n")
    
    o.write(f"""python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/add_outgroup_to_vcf.py D_SIMULANS D_MELANOGASTER \\\n""")
    o.write(f"""/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_mel.maf \\\n""")
    o.write(f"""/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/mel_aligned_sim.fasta \\\n""")
    o.write(f"""/work/users/d/t/dturissi/drosophila/ssh/genotyping/filtered_vcfs/sim_sech_hybrid_filtered_biallelic_snp_${{chroms[${{i}}]}}.vcf.gz \\\n""")
    o.write(f"""/work/users/d/t/dturissi/drosophila/ssh/genotyping/filtered_vcfs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup_${{chroms[${{i}}]}}.vcf.gz\n""")
      
       
if __name__ == '__main__':
  main()
