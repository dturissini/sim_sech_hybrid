#!/bin/bash
#SBATCH -p general
#sbatch -J canu_sech
#SBATCH --cpus-per-task=1
#SBATCH -t 144:00:00
#SBATCH --mem=8g
#SBATCH --array=1-2
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/canu_sech_genomes_low_cov_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/canu_sech_genomes_low_cov_%j.err


sech_sample_ids=('sech_Denis_NF45_female' 'sech_Denis_papaya152')

i=$((${SLURM_ARRAY_TASK_ID} - 1))
python3 assemble_sech_genome_low_cov.py ${sech_sample_ids[${i}]}

