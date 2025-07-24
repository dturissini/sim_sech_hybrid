import os
import sys


#python3 make_biallelic_snp_vcf_to_zarr_script.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  log_dir = os.path.join(base_dir, 'logs')  
  filtered_vcf_dir = os.path.join(base_dir, 'filtered_vcfs')
  zarr_dir = os.path.join(base_dir, 'zarrs')
    
  vcf_to_zarr_sbatch_file = os.path.join(base_dir, 'chrom_biallelic_snp_vcf_to_zarr.sbatch')
  
  os.system(f"mkdir -p {zarr_dir}")

  
  chroms = ['2L', '2R', '3L', '3R', 'X', '4']
                                                                            
  total_cmds = 0
  chrom_str = ''
  vcf_str = ''                            
  zarr_str = ''  
  for chrom in chroms:
    total_cmds += 1

    chrom_filtered_vcf_file = os.path.join(filtered_vcf_dir, 'sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup_' + chrom + '.vcf.gz')          
    chrom_zarr_file = os.path.join(zarr_dir, 'sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup_' + chrom + '.zarr')
        
    chrom_str += '"' + chrom + '" '
    vcf_str += '"' + chrom_filtered_vcf_file + '" '
    zarr_str += '"' + chrom_zarr_file + '" '


  chrom_str = chrom_str.rstrip()
  vcf_str = vcf_str.rstrip()
  zarr_str = zarr_str.rstrip()

  with open(vcf_to_zarr_sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J vcf_to_zarr\n")
    o.write(f"#SBATCH -t 24:00:00\n")
    o.write(f"#SBATCH --mem=16g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}%100\n")
    o.write(f"#SBATCH -o {log_dir}/biallelic_snp_vcf_to_zarr_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/biallelic_snp_vcf_to_zarr_%A_%a.err\n\n")
    
    
    o.write(f"""ml python/3.12.2\n\n\n""")
    
    o.write(f"""chroms=({chrom_str})\n\n""")
    o.write(f"""vcfs=({vcf_str})\n\n""")
    o.write(f"""zarrs=({zarr_str})\n\n""")

    o.write(f"i=$((${{SLURM_ARRAY_TASK_ID}} - 1))\n")
    
    o.write(f"""python3 {base_dir}/filtered_vcf_to_zarr.py ${{chroms[${{i}}]}} ${{vcfs[${{i}}]}} ${{zarrs[${{i}}]}}\n""")


       
if __name__ == '__main__':
  main()
