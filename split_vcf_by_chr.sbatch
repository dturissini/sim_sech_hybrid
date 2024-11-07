#!/bin/bash
#SBATCH -p general
#sbatch -J split_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=4g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/split_vcf_by_chr_%j_%a.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/split_vcf_by_chr_%j_%a.err

chroms=('2L' '2R' '3L' '3R' 'X')

i=$((${SLURM_ARRAY_TASK_ID} - 1))

vcftools --recode --maf .01 --gzvcf /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz --chr ${chroms[${i}]} --out /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_${chroms[${i}]}
gzip /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_${chroms[${i}]}.recode.vcf
