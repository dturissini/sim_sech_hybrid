#!/bin/bash
#SBATCH -p general
#sbatch -J hal2maf
#SBATCH --cpus-per-task=1
#SBATCH -t 8:00:00
#SBATCH --mem=8g
#SBATCH -o /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_freq_dist_%j.out
#SBATCH -e /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_freq_dist_%j.err


source /proj/matutelb/software/cactus-bin-v2.8.0/venv-cactus-v2.8.0/bin/activate
python3 /proj/matutelb/projects/gwas/scripts/get_sim_sech_freq_dist.py /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.vcf.gz /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_final_samplelist.tsv /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_freq_dist.db
