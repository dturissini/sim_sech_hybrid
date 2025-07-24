import os
import sys


#python3 make_filtered_vcf_by_chrom_script.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  vcf_dir = os.path.join(base_dir, 'vcfs')
  log_dir = os.path.join(base_dir, 'logs')  
  filtered_vcf_dir = os.path.join(base_dir, 'filtered_vcfs')
  
  chr_mapping_file = os.path.join(base_dir, 'chr_mapping_file.txt')
  sbatch_file = os.path.join(base_dir, 'filtered_vcf_by_chrom.sbatch')
  chrom_rename_file = os.path.join(base_dir, 'simulans_chrom_rename.txt')

  filtered_vcf_file = os.path.join(filtered_vcf_dir, 'sim_sech_hybrid_filtered.vcf.gz')
  
  lk_chrom = {'NC_052520.2': '2L',
              'NC_052521.2': '2R',
              'NC_052522.2': '3L',
              'NC_052523.2': '3R',
              'NC_052525.2': 'X',
              'NC_052524.2': '4'}
              
           
  total_cmds = len(lk_chrom)

  chrom_str = ''
  raw_chrom_str = ''
  with open(chrom_rename_file, 'w') as r: 
    for raw_chrom in lk_chrom:
      r.write(f"{raw_chrom} {lk_chrom[raw_chrom]}\n")      

      chrom_str += '"' + lk_chrom[raw_chrom] + '" '
      raw_chrom_str += '"' + raw_chrom + '" '
 

  chrom_str = chrom_str.rstrip()
  raw_chrom_str = raw_chrom_str.rstrip()

                              
  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J filter_chrom_vcf\n")
    o.write(f"#SBATCH -t 96:00:00\n")
    o.write(f"#SBATCH --mem=8g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}\n")
    o.write(f"#SBATCH -o {log_dir}/drosophila_filter_vcf_by_chrom_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/drosophila_filter_vcf_by_chrom_%A_%a.err\n\n")
    
    o.write(f"ml samtools/1.21\n")

    o.write(f"""chroms=({chrom_str})\n\n""")
    o.write(f"""raw_chroms=({raw_chrom_str})\n\n""")

    o.write(f"i=$((${{SLURM_ARRAY_TASK_ID}} - 1))\n\n")
    
    o.write(f"""chrom_vcf_file_old_name=/work/users/d/t/dturissi/drosophila/ssh/genotyping/filtered_vcfs/sim_sech_hybrid_filtered_biallelic_snp_${{chroms[${{i}}]}}_old_chrom_name.vcf.gz\n""")
    o.write(f"""chrom_vcf_file=/work/users/d/t/dturissi/drosophila/ssh/genotyping/filtered_vcfs/sim_sech_hybrid_filtered_biallelic_snp_${{chroms[${{i}}]}}.vcf.gz\n\n""")
    
    o.write(f"""bcftools view {filtered_vcf_file} -V indels -f .,PASS -M 2 -r ${{raw_chroms[${{i}}]}} -O v -o ${{chrom_vcf_file_old_name}}\n""")
    o.write(f"""bcftools annotate ${{chrom_vcf_file_old_name}} --rename-chrs {chrom_rename_file} -O v -o ${{chrom_vcf_file}}\n""")
      

       
if __name__ == '__main__':
  main()
