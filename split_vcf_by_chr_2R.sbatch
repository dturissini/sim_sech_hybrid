#!/bin/bash
#SBATCH -p general
#sbatch -J split_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 4:00:00
#SBATCH --mem=4g
#SBATCH -o /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/split_vcf_by_chr_2R_%j.out
#SBATCH -e /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/split_vcf_by_chr_2R_%j.err

vcftools --recode --maf .01 --gzvcf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz --chr 2R --out /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2R

gzip /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2R.recode.vcf
