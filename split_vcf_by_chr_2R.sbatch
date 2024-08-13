#!/bin/bash
#SBATCH -p general
#sbatch -J split_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 4:00:00
#SBATCH --mem=4g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/split_vcf_by_chr_2R_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/split_vcf_by_chr_2R_%j.err

vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr 2R --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2R

gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2R.recode.vcf
