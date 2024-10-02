import os
import sys
import tempfile

#python3 inversion_check_lastz_reads.py sech_Denis_NF69



def main():
  sample_id = sys.argv[1]
  
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly'
  assembly_dir = os.path.join(base_dir, 'assemblies', sample_id)
  inv_check_dir = os.path.join(base_dir, 'inversion_check')
  inv_overlap_dir =  os.path.join(inv_check_dir, 'inversion_overlaps_trimmedReads')
  maf_dir = os.path.join(inv_check_dir, 'mafs_reads')

  fasta_file = os.path.join(assembly_dir, sample_id + '.trimmedReads.fasta')
  inv_overlap_file = os.path.join(inv_overlap_dir, sample_id + '_inv_overlap_trimmedReads.txt')

  ref_genome = '/work/users/d/t/dturissi/drosophila/ssh/assembly/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'
  inversion_file = os.path.join(inv_check_dir, 'denis_possible_inversions.txt')
  
  os.system(f"mkdir -p {maf_dir}")
  
  #unzip trimmedReds file if zipped since lastz can't handle gz files
  if os.path.isfile(f"{fasta_file}.gz"):
    os.system(f"gunzip {fasta_file}.gz")
  
  #process and store inversion locations
  inversions = {}
  with open(inversion_file, 'r') as i:
    next(i)
    for line in i:
      inv_name, chrom, start, end = line.strip().split()
      inversions[inv_name] = {'chrom': chrom, 'start': start, 'end': end}
      
  with open(inv_overlap_file, 'r') as io:
    for line in io:
      inv_name, read, read_len, read_start, read_end, chrom, chrom_start, chrom_end, strand = line.strip().split()    
      
      if (int(chrom_end) - int(chrom_start) + 1 > 1000):
        print(inv_name, read)
          
        with tempfile.NamedTemporaryFile(mode='w') as target_chrom:
          with tempfile.NamedTemporaryFile(mode='w') as query_read: 
            target_chrom.write(f"{chrom}")
            query_read.write(f"{read}")
            
            target_chrom.flush() 
            query_read.flush() 
            
            #run lastz
            maf_file = os.path.join(maf_dir, inv_name + '_' + sample_id + '_' + read + '.maf')  
            os.system(f"""lastz {ref_genome}[subset={target_chrom.name}][{inversions[inv_name]['start']}..{inversions[inv_name]['end']}] {fasta_file}[subset={query_read.name}] --notransition --step=20 --nogapped --format=maf > {maf_file}""")
  


if __name__ == '__main__':
  main()
