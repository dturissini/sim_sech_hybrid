import sqlite3
import os
import sys


#python3 make_map_reads_scripts.py



def main():  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/genotyping'
  fastq_dir = os.path.join(base_dir, 'fastqs')
  map_dir = os.path.join(base_dir, 'mapping')
  map_sh_dir = os.path.join(map_dir, 'map_sh')
  bam_dir = os.path.join(map_dir, 'bams')
  log_dir = os.path.join(map_dir, 'logs')

  ref_genome_file = os.path.join(base_dir, 'ref_genome', 'GCF_016746395.2_Prin_Dsim_3.1_genomic.fna')
  db_file = os.path.join(base_dir, 'ssh_genotyping.db')
  
  
  threads = 8

  master_sh_file = os.path.join(map_dir, 'map_ssh_reads.sh')
  sbatch_file = os.path.join(map_dir, 'map_ssh_reads.sbatch')
 
  os.system(f"mkdir -p {map_sh_dir}")
  os.system(f"mkdir -p {bam_dir}")
  os.system(f"mkdir -p {log_dir}")
  
  conn = sqlite3.connect(db_file)    
  sample_query = conn.execute(f"""select distinct sample_id
                                  from ssh_fastqs
                                  order by sample_id""")
                              
  total_cmds = 0
  with open(master_sh_file, 'w') as m:
    for sample_id, in sample_query:
      total_cmds += 1
      
      bam_sample_dir = os.path.join(bam_dir, sample_id)
      sample_sh_file = os.path.join(map_sh_dir, sample_id + '_map_reads.sh')
      bam_file_base = os.path.join(bam_sample_dir, sample_id)

      os.system(f"mkdir -p {bam_sample_dir}")
      
      m.write(f"{sample_sh_file}\n")

      sample_bams = []
      with open(sample_sh_file, 'w') as i:
        i.write(f"set -o pipefail\n\n")
        
        fastq_base_query = conn.execute(f"""select fastq_base
                                            from ssh_fastqs
                                            where sample_id = '{sample_id}'
                                            group by fastq_base""")
      
        for fastq_base, in fastq_base_query:
          bam_file = fastq_base + '.bam'
          samtools_sort_tmp_file = os.path.join(bam_sample_dir, fastq_base + '_samtools_tmp')
          sample_bams.append(bam_file)
          
          rg_str = '@RG\\tID:' + sample_id + '\\tSM:' + sample_id + '\\tLB:ga\\tPL:Illumina'
          
          fastq_query = conn.execute(f"""select fastq_file
                                         from ssh_fastqs f
                                         where sample_id = '{sample_id}'
                                         and fastq_base = '{fastq_base}'
                                         group by fastq_file""")
          
          fastq_files = []
          for fastq_file, in fastq_query:
            fastq_files.append(fastq_file) 
          
          if len(fastq_files) == 1:
            i.write(f"""bwa-mem2 mem -t {threads} -Y -R "{rg_str}" {ref_genome_file} {fastq_files[0]}  | \\\n""")
            i.write(f"""samtools fixmate -u -m - - | \\\n""")
            i.write(f"""samtools sort -u -@{threads} -T {samtools_sort_tmp_file} - | \\\n""")
            i.write(f"""samtools markdup -@{threads} --reference {ref_genome_file} -O bam - {bam_file}\n\n""")
          else:
            i.write(f"""bwa-mem2 mem -t {threads} -Y -R "{rg_str}" {ref_genome_file} {fastq_files[0]} {fastq_files[1]} | \\\n""")
            i.write(f"""samtools fixmate -u -m - - | \\\n""")
            i.write(f"""samtools sort -u -@{threads} -T {samtools_sort_tmp_file} - | \\\n""")
            i.write(f"""samtools markdup -@{threads} --reference {ref_genome_file} -O bam - {bam_file}\n\n""")
            

#rename and merge bams
        merged_bam_file = bam_file_base + '.bam'
        if len(sample_bams) == 1:
          i.write(f"mv {sample_bams[0]} {merged_bam_file}\n\n")
        else:
          sample_bams_str = ' '.join(sample_bams)
          i.write(f"samtools merge {merged_bam_file} {sample_bams_str}\n\n")
          i.write(f"rm -r {sample_bams_str}\n\n")

        i.write(f"samtools index {merged_bam_file}\n\n")

      os.system(f"chmod +x {sample_sh_file}")



  with open(sbatch_file, 'w') as o:
    o.write(f"#!/bin/bash\n")
    o.write(f"#SBATCH -p general\n")
    o.write(f"#sbatch -J ssh_bwa\n")
    o.write(f"#SBATCH --cpus-per-task={threads}\n")
    o.write(f"#SBATCH -t 168:00:00\n")
    o.write(f"#SBATCH --mem=32g\n")
    o.write(f"#SBATCH --array=1-{total_cmds}%50\n")
    o.write(f"#SBATCH -o {log_dir}/ssh_map_reads_%A_%a.out\n")
    o.write(f"#SBATCH -e {log_dir}/ssh_map_reads_%A_%a.err\n\n")
    
    o.write(f"ml bwa-mem2/2.2.1\n")
    o.write(f"ml samtools/1.21\n\n")

    o.write(f"bwa_sh=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p {master_sh_file})\n")
    o.write(f"""eval "${{bwa_sh}}"\n""")


       
if __name__ == '__main__':
  main()
