import argparse
import os

#python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/make_outgroup_fasta_from_maf.py D_SIMULANS D_MELANOGASTER /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/sim_mel.maf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/mel_aligned_sim.fasta 


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_species")
    parser.add_argument("outgroup_species")
    parser.add_argument("maf_file")
    parser.add_argument("fasta_file")

    return parser.parse_args()


def main():
  args = parse_args()
  ref_species = args.ref_species  
  outgroup_species = args.outgroup_species  
  maf_file = args.maf_file
  fasta_file = args.fasta_file


  min_depth = 6
  lk_chrom = {'NC_052520.2': '2L',
              'NC_052521.2': '2R',
              'NC_052522.2': '3L',
              'NC_052523.2': '3R',
              'NC_052524.2': '4',
              'NC_052525.2': 'X'}
  
  true_len = {'2L':23857595,
              '2R':22319025,
              '3L':23399903,
              '3R':28149585,
              '4' :1146867,
              'X' :22032822}
    
   
  outgroup_syn_chroms = {}            
  
  block_ref_start = 0
  block_ref_end = 0
  block_ref_chrom = ''
  block_ref_insertions = []
  last_chrom = ''
  counter = 0
  with open(maf_file, 'r') as m: 
    for line in m:
      if counter % 100000 == 0:
        print(counter)
      counter += 1

      line = line.strip()
      if line[:1] != '#':
        if line[:1] in ['s']:
          (dummy, src, start, size, strand, chrom_len, seq) = line.split('\t')
          
          (species, chrom) = src.split('.')
                    
          if species == ref_species:
            if chrom not in outgroup_syn_chroms:
              outgroup_syn_chroms[chrom] = 'N' * int(chrom_len)
            
            if chrom != last_chrom:
              print(chrom, last_chrom, outgroup_syn_chroms.keys())
            
            block_ref_insert_idx = [i for i, nt in enumerate(seq) if nt == '-']
            seq = seq.replace('-', '')
            
            block_ref_start = int(start)
            block_ref_end = block_ref_start + len(seq)
            block_ref_chrom = chrom  
            last_chrom = chrom
          elif species == outgroup_species:
            if block_ref_insert_idx:
              seq = "".join([nt for idx, nt in enumerate(seq) if idx not in block_ref_insert_idx])
            
            outgroup_syn_chroms[block_ref_chrom] = outgroup_syn_chroms[block_ref_chrom][:block_ref_start] + seq.upper() + outgroup_syn_chroms[block_ref_chrom][block_ref_end:]
            
            
            if block_ref_end - block_ref_start != len(seq):
              print(counter, 'seq wrong len', block_ref_chrom, block_ref_start, block_ref_end, len(seq), line)
            
            if block_ref_chrom in true_len:
              if true_len[block_ref_chrom] != len(outgroup_syn_chroms[block_ref_chrom]):
                print(counter, 'chroms diff', block_ref_chrom, true_len[block_ref_chrom], len(outgroup_syn_chroms[block_ref_chrom]), line)
            

  with open(fasta_file, 'w') as o: 
    for chrom in sorted(outgroup_syn_chroms.keys()):
      print(chrom, len(outgroup_syn_chroms[chrom]))
      o.write(f">{chrom}\n{outgroup_syn_chroms[chrom]}\n")

    

if __name__ == '__main__':
  main()
