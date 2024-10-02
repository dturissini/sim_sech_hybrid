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
  inv_overlap_dir = os.path.join(inv_check_dir, 'inversion_overlaps_' + seq_type)

  gz_suffix = ''
  if seq_type == 'trimmedReads':
    gz_suffix = '.gz'
   
  fasta_file = os.path.join(assembly_dir, sample_id + '.' + seq_type + '.fasta' + gz_suffix)
  paf_file = os.path.join(paf_dir, sample_id + '.paf')
  sam_file = os.path.join(sam_dir, sample_id + '.sam')
  inv_overlap_file = os.path.join(inv_overlap_dir, sample_id + '_inv_overlap_' + seq_type + '.txt')
  
  ref_genome = '/work/users/d/t/dturissi/drosophila/ssh/assembly/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna'
  
  os.system(f"mkdir -p {paf_dir}")
  os.system(f"mkdir -p {sam_dir}")
  os.system(f"mkdir -p {inv_overlap_dir}")
  
  
  #run minimap
  os.system(f"minimap2 -x asm20 {ref_genome} {fasta_file} > {paf_file}")
  os.system(f"minimap2 -x asm20 -a {ref_genome} {fasta_file} > {sam_file}")
  
  
  inversion_file = 'denis_possible_inversions.txt'
  inversions = {}
  with open(inversion_file, 'r') as i:
    next(i)
    for line in i:
      (name, chrom, start, end) = line.split("\t")
      inversions[name] = [chrom, int(start), int(end)]
  

  #process putative inversions
  counter = 0
  with open(inv_overlap_file, 'w') as o:
    with open(paf_file, 'r') as p:
      for line in p:
        if counter % 100000 == 0:
          {print(counter)}
        counter += 1
        
        values = line.split("\t")
        contig_name = values[0]
        contig_len = values[1]
        contig_start = int(values[2]) + 1
        contig_end = int(values[3]) + 1
        strand = values[4]
        chrom = values[5]
        start = int(values[7]) + 1
        end = int(values[8]) + 1
        
        for inv_name in inversions:
          if inversions[inv_name][0] == chrom:
            if (inversions[inv_name][1] >= start and inversions[inv_name][1] <= end) or (inversions[inv_name][2] >= start and inversions[inv_name][2] <= end) or (start >= inversions[inv_name][1] and start <= inversions[inv_name][2]) or (end >= inversions[inv_name][1] and end <= inversions[inv_name][2]):
              o.write(f"{inv_name}\t{contig_name}\t{contig_len}\t{contig_start}\t{contig_end}\t{chrom}\t{start}\t{end}\t{strand}\n")



if __name__ == '__main__':
  main()
