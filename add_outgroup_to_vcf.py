import argparse
import os
import re
import gzip

#python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/add_outgroup_to_vcf.py D_SIMULANS D_MELANOGASTER /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_mel.maf /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/mel_aligned_sim.fasta /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.vcf.gz /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_species")
    parser.add_argument("outgroup_species")
    parser.add_argument("maf_file")
    parser.add_argument("outgroup_fasta_file")    
    parser.add_argument("vcf_file_in")
    parser.add_argument("vcf_file_out")

    return parser.parse_args()


def main():
  args = parse_args()
  ref_species = args.ref_species  
  outgroup_species = args.outgroup_species  
  maf_file = args.maf_file
  outgroup_fasta_file = args.outgroup_fasta_file
  vcf_file_in = args.vcf_file_in
  vcf_file_out = args.vcf_file_out  

  vcf_file_out_plink = vcf_file_out[:-7] + '_plink.vcf.gz'
  vcf_file_out_multi = vcf_file_out[:-7] + '_multi.vcf.gz'
  vcf_file_out_no_mel = vcf_file_out[:-7] + '_no_mel.vcf.gz'

  min_depth = 6
  lk_chrom = {'NC_052520.2': '2L',
              'NC_052521.2': '2R',
              'NC_052522.2': '3L',
              'NC_052523.2': '3R',
              'NC_052525.2': 'X',
              'NC_052524.2': '4'}
  
  chrom_plink_rename = {'NC_052520.2': '1',
                        'NC_052521.2': '2',
                        'NC_052522.2': '3',
                        'NC_052523.2': '4',
                        'NC_052525.2': '5',
                        'NC_052524.2': '6'}
  
 
  outgroup_syn_chroms = {}
  chrom = ''            
  with open(outgroup_fasta_file, 'r') as f:
    for line in f:
      line = line.strip()
      if line[:1] == '>':
        chrom = line[1:]
      else: 
        outgroup_syn_chroms[chrom] = line
   
  counter = 0
  with gzip.open(vcf_file_out_plink, 'wt') as op: 
    with gzip.open(vcf_file_out_multi, 'wt') as om: 
      with gzip.open(vcf_file_out_no_mel, 'wt') as onm: 
        with gzip.open(vcf_file_out, 'wt') as o: 
          with gzip.open(vcf_file_in, 'rt') as i: 
            for line in i:
              if counter % 100000 == 0:
                print(counter)
              
              counter += 1
        
              if line[:2] == "##":
                o.write(line) 
                op.write(line) 
                om.write(line) 
                onm.write(line) 
              elif line[:1] == "#":
                header = line.split('\t')   
                header = header[:9] + [outgroup_species] + header[9:]      
                header_out = "\t".join(header)     
                o.write(f"{header_out}")     
                op.write(f"{header_out}")     
                om.write(f"{header_out}")     
                onm.write(f"{header}")     
              else:
                values = line.split('\t')       
                if values[0] in lk_chrom:   
                  sim_chrom = lk_chrom[values[0]]
                  plink_chrom = chrom_plink_rename[values[0]]
                  mel_allele = outgroup_syn_chroms[sim_chrom][(int(values[1]) - 1):int(values[1])]
                  if mel_allele not in ['N', '-']:
                    alleles = {values[3], values[4], mel_allele}
                    if len(alleles) == 2:
                      outgroup_GT = "0/0"
                      if mel_allele == values[4]:
                        outgroup_GT = "1/1"
                        
                      outgroup_genotype = outgroup_GT + ":99:99"
                      values = values[:9] + [outgroup_genotype] + values[9:]       
                      line_out = "\t".join(values)     

                      values[0] = sim_chrom
                      line_out = "\t".join(values)     
                      o.write(f"{line_out}")   

                      values[0] = plink_chrom
                      line_out = "\t".join(values)     
                      op.write(f"{line_out}")   
                    else:
                      #eventually add support for mel introducing third (or fourth) allele
                      om.write(f"{line}")   
                  else:
                     onm.write(f"{line}")   
    

if __name__ == '__main__':
  main()
