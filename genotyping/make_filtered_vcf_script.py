import os
import sys


#python3 make_filtered_vcf_script.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'

  vcf_dir = os.path.join(base_dir, 'vcfs')
  filtered_vcf_dir = os.path.join(base_dir, 'filtered_vcfs')
  log_dir = os.path.join(base_dir, 'logs')
  
  vcf_file = os.path.join(vcf_dir, 'sim_sech_hybrid.vcf.gz')
  filtered_vcf_file = os.path.join(filtered_vcf_dir, 'sim_sech_hybrid_filtered.vcf.gz')

  ref_genome_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'    
  sbatch_file = os.path.join(base_dir, 'make_filtered_vcf.sbatch')
 
  os.system(f"mkdir -p {filtered_vcf_dir}")

  
                              
  with open(sbatch_file, 'w') as sm:    
    sm.write(f"#!/bin/bash\n")
    sm.write(f"#SBATCH -p general\n")
    sm.write(f"#sbatch -J filtered_vcf\n")
    sm.write(f"#SBATCH -t 96:00:00\n")
    sm.write(f"#SBATCH --mem=4g\n")
    sm.write(f"#SBATCH -o {log_dir}/filtered_vcf_%A_%a.out\n")
    sm.write(f"#SBATCH -e {log_dir}/filtered_vcf_%A_%a.err\n\n")
    
    sm.write(f"ml gatk/4.6.1.0\n")

    sm.write(f"""gatk VariantFiltration \\\n""")
    sm.write(f"""-R {ref_genome_file} \\\n""")
    sm.write(f"""-V {vcf_file} \\\n""")
    sm.write(f"""-O {filtered_vcf_file} \\\n""")
    sm.write(f"""-filter-name "QD_filter" -filter "QD < 2.0" \\\n""")
    sm.write(f"""-filter-name "FS_filter" -filter "FS > 60.0" \\\n""")
    sm.write(f"""-filter-name "MQ_filter" -filter "MQ < 40.0" \\\n""")
    sm.write(f"""-filter-name "SOR_filter" -filter "SOR > 4.0" \\\n""")
    sm.write(f"""-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \\\n""")
    sm.write(f"""-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"\n\n""")
                            
            

       
if __name__ == '__main__':
  main()
