#!/bin/bash
#SBATCH -p general
#sbatch -J split_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=4g
#SBATCH -o /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/split_vcf_by_chr_%j.out
#SBATCH -e /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/split_vcf_by_chr_%j.err

vcftools --recode --maf .01 --gzvcf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz --chr 2L --out /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2L
vcftools --recode --maf .01 --gzvcf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz --chr 2R --out /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2R
vcftools --recode --maf .01 --gzvcf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz --chr 3L --out /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3L
vcftools --recode --maf .01 --gzvcf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz --chr 3R --out /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3R
vcftools --recode --maf .01 --gzvcf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz --chr X  --out /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_X 

gzip /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2L.recode.vcf
gzip /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2R.recode.vcf
gzip /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3L.recode.vcf
gzip /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3R.recode.vcf
gzip /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_X.recode.vcf
