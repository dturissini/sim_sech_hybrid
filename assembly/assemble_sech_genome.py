import os
import tempfile
import sys


#python3 assemble_sech_genome.py sech_Denis_NF22



def main():
  sech_sample_id = sys.argv[1]
  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly'
  fastq_dir = os.path.join(base_dir, 'fastqs')
  fastq_file = os.path.join(fastq_dir, sech_sample_id + '.fastq')
  original_fastq_file = os.path.join('/proj/matutelb/projects/sech_7_19_24', sech_sample_id + '.fastq.gz')
  
  assembly_dir = os.path.join(base_dir, 'assemblies', sech_sample_id)
  
  
  #mkdir fastqs
  if 'female' in sech_sample_id:
    os.system(f"gunzip -c {original_fastq_file} > {fastq_file}")
  else:
    os.system(f"cp {original_fastq_file} {fastq_file}")
  

  #make sample dir
  os.system(f"mkdir -p {assembly_dir}")

  os.system(f"""/proj/matutelb/software/canu-2.2/bin/canu -p {sech_sample_id} \
                                          -d {assembly_dir} \
                                          genomeSize=130m \
                                          -nanopore {fastq_file} \
                                          -gridOptions="--time=144:00:00 --mem-per-cpu=4g --partition=general" """)


if __name__ == '__main__':
  main()
