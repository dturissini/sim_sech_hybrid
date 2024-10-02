#!/bin/bash
#SBATCH -p general
#sbatch -J canu_sech
#SBATCH --cpus-per-task=1
#SBATCH -t 144:00:00
#SBATCH --mem=8g
#SBATCH --array=1-1
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/canu_sech_genomes_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/canu_sech_genomes_%j.err


sech_sample_ids=('sech_Denis_NF69')
##sech_sample_ids=('sech_Denis_NF69' 'sech_Anro_114_female')

##sech_sample_ids=('sech_Denis_NF22' 'sech_Denis_NF69' 'sech_Denis_NF712' 'sech_Denis_papaya152' 'sech_Anro_114_female' 'sech_Denis_NF45_female')

i=$((${SLURM_ARRAY_TASK_ID} - 1))
python3 assemble_sech_genome.py ${sech_sample_ids[${i}]}

