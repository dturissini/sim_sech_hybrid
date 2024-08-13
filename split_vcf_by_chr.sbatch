#!/bin/bash
#SBATCH -p general
#sbatch -J split_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=4g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/split_vcf_by_chr_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/split_vcf_by_chr_%j.err

vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr 2L --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2L
vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr 2R --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2R
vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr 3L --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3L
vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr 3R --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3R
vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr X  --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_X 

gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2L.recode.vcf
gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2R.recode.vcf
gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3L.recode.vcf
gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3R.recode.vcf
gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_X.recode.vcf
