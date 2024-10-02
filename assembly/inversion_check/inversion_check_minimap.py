import os
import tempfile
import sys


#python3 inversion_check_minimap.py sech_Denis_NF49 contigs
#python3 inversion_check_minimap.py sech_Denis_NF49 trimmedReads



def main():
  sample_id = sys.argv[1]
  seq_type = sys.argv[2]
  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly'
  assembly_dir = os.path.join(base_dir, 'assemblies', sample_id)
  inv_check_dir = os.path.join(base_dir, 'inversion_check')
  paf_dir = os.path.join(inv_check_dir, 'pafs')
  sam_dir = os.path.join(inv_check_dir, 'sams')

  gz_suffix = ''
  if seq_type == 'trimmedReads':
    gz_suffix = '.gz'
   
  fasta_file = os.path.join(assembly_dir, sample_id + '.' + seq_type + '.fasta' + gz_suffix)
  paf_file = os.path.join(paf_dir, sample_id + '.paf')
  sam_file = os.path.join(sam_dir, sample_id + '.sam')
  
  ref_genome = '/work/users/d/t/dturissi/drosophila/ssh/assembly/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'
  
  os.system(f"mkdir -p {paf_dir}")
  os.system(f"mkdir -p {sam_dir}")
  
  
  #run minimap
  os.system(f"minimap2 -x asm20 {ref_genome} {fasta_file} > {paf_file}")
  os.system(f"minimap2 -x asm20 -a {ref_genome} {fasta_file} > {sam_file}")
  
  

if __name__ == '__main__':
  main()
