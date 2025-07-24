import sqlite3
import os
import sys


#python3 make_vcf_scripts.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  
  gvcf_dir = os.path.join(base_dir, 'gvcfs')
  log_dir = os.path.join(base_dir, 'logs')
  genotype_tmp_dir = os.path.join(base_dir, 'tmp_gatk_GenotypeGVCFs')
  db_dir = os.path.join(base_dir, 'dbs')
  vcf_dir = os.path.join(base_dir, 'vcfs')
  sample_map_dir = os.path.join(base_dir, 'sample_maps')
  
  gatkdb_file = os.path.join(db_dir, 'sim_sech_hybrid.gatkdb')
  sample_map_file = os.path.join(sample_map_dir, 'sim_sech_hybrid.sample_map')
  vcf_file = os.path.join(vcf_dir, 'sim_sech_hybrid.vcf.gz')
  gvcf_file = os.path.join(gvcf_dir, 'sim_sech_hybrid.gvcf.gz')      
  
  
  ref_genome_file = '/work/users/d/t/dturissi/drosophila/ssh/genotyping/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'  
  interval_list_file = ref_genome_file + '.list'
  
  
  db_sbatch_file = os.path.join(base_dir, 'make_gatk_dbs.sbatch')
  vcf_sbatch_file = os.path.join(base_dir, 'make_vcfs.sbatch')
 
  os.system(f"mkdir -p {vcf_dir}")
  os.system(f"mkdir -p {gvcf_dir}")
  os.system(f"mkdir -p {db_dir}")
  os.system(f"mkdir -p {sample_map_dir}")
  os.system(f"mkdir -p {log_dir}")
  os.system(f"mkdir -p {genotype_tmp_dir}")
  
  threads = 4
  
  db_file = os.path.join(base_dir, 'ssh_genotyping.db')
  conn = sqlite3.connect(db_file)    
  sample_query = conn.execute(f"""select sample_id
                                  from bam_files
                                  order by sample_id""")


  with open(sample_map_file, 'w') as s:
    for sample_id, in sample_query:
      gvcf_sample_dir = os.path.join(gvcf_dir, sample_id)
      gvcf_file = os.path.join(gvcf_sample_dir, sample_id + '.gvcf.gz')
      gvcf_tbi_file = gvcf_file + '.tbi'
      
      s.write(f"{sample_id}\t{gvcf_file}\t{gvcf_tbi_file}\n")   
                              
                                                  
    

  with open(db_sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J ssh_gatk_db\n")
    o.write(f"#SBATCH -t 240:00:00\n")
    o.write(f"#SBATCH --mem=64g\n")
    o.write(f"#SBATCH --cpus-per-task={threads}\n")
    o.write(f"#SBATCH -o {log_dir}/ssh_gatk_db_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/ssh_gatk_db_%A_%a.err\n\n")
    
    o.write(f"ml gatk/4.6.1.0\n")

    o.write(f"""rm -rf {gatkdb_file}\n\n""")
    o.write(f"""gatk --java-options "-Xmx60g -Xms60g" \\\n""")
    o.write(f"""GenomicsDBImport \\\n""")
    o.write(f"""--genomicsdb-workspace-path {gatkdb_file} \\\n""")
    o.write(f"""--merge-contigs-into-num-partitions 10 \\\n""")
    o.write(f"""--batch-size 20 \\\n""")
    o.write(f"""-L {interval_list_file} \\\n""")
    o.write(f"""--sample-name-map {sample_map_file} \\\n""")
    o.write(f"""--reader-threads {threads}\n\n""")    


  with open(vcf_sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J ssh_vcf\n")
    o.write(f"#SBATCH -t 240:00:00\n")
    o.write(f"#SBATCH --mem=32g\n")
    o.write(f"#SBATCH -o {log_dir}/ssh_vcf_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/ssh_vcf_%A_%a.err\n\n")
    
    o.write(f"ml gatk/4.6.1.0\n")

    o.write(f"""gatk --java-options "-Xmx30g" SelectVariants \\\n""")
    o.write(f"""-R {ref_genome_file} \\\n""")
    o.write(f"""-V gendb://{gatkdb_file} \\\n""")
    o.write(f"""-O {gvcf_file} \\\n""")
    o.write(f"""--tmp-dir {genotype_tmp_dir}\n""")
    
    o.write(f"""\n\n""")
    
    o.write(f"""gatk --java-options "-Xmx30g" GenotypeGVCFs \\\n""")
    o.write(f"""-R {ref_genome_file} \\\n""")
    o.write(f"""-V {gvcf_file} \\\n""")
    o.write(f"""-O {vcf_file} \\\n""")
    o.write(f"""-all-sites \\\n""")
    o.write(f"""--tmp-dir {genotype_tmp_dir}\n""")




#old code to directly call all-sites genotype from GenomicsDBImport
#the direct approach is untenable due to a known memory leak in gatk:
#https://github.com/broadinstitute/gatk/issues/8989

##  with open(vcf_sbatch_file, 'w') as o:
##    o.write(f"#!/bin/bash\n")
##    o.write(f"#SBATCH -p general\n")
##    o.write(f"#sbatch -J vcf\n")
##    o.write(f"#SBATCH -t 240:00:00\n")
##    o.write(f"#SBATCH --mem=154g\n")
##    o.write(f"#SBATCH --array=1-{total_cmds}%100\n")
##    o.write(f"#SBATCH -o {log_dir}/drosophila_raw_vcf_%A_%a.out\n")
##    o.write(f"#SBATCH -e {log_dir}/drosophila_raw_vcf_%A_%a.err\n\n")
##    
##    o.write(f"ml gatk/4.6.1.0\n")
##
##    o.write(f"gatk_sh=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {vcf_master_sh_file})\n")
##    o.write(f"""eval "${{gatk_sh}}"\n""")


 

#old code to directly call all-sites genotype from GenomicsDBImport
#the direct approach is untenable due to a known memory leak in gatk:
#https://github.com/broadinstitute/gatk/issues/8989
##        with open(vcf_sh_file, 'w') as i:    
##          i.write(f"""gatk --java-options "-Xmx150g" GenotypeGVCFs \\\n""")
##          i.write(f"""-R {ref_genome_file} \\\n""")
##          i.write(f"""-V gendb://{gatkdb_file} \\\n""")
##          i.write(f"""-O {vcf_file} \\\n""")
##          i.write(f"""-all-sites \\\n""")
##          i.write(f"""--tmp-dir {genotype_tmp_dir}\n""")

            


       
if __name__ == '__main__':
  main()
