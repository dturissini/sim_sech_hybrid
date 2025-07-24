#!/bin/bash
#SBATCH -p general
#sbatch -J canu_sech
#SBATCH --cpus-per-task=1
#SBATCH -t 144:00:00
#SBATCH --mem=8g
#SBATCH --array=1-3
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/canu_sech_genomes_%A_%a.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/canu_sech_genomes_%A_%a.err


##sech_sample_ids=('sech_Denis_NF22' 'sech_Denis_NF69' 'sech_Denis_NF712' 'sech_Denis_papaya152' 'sech_Anro_114_female' 'sech_Denis_NF45_female')
sech_sample_ids=('DenisNF134' 'SECH_20-sech_PNF5_GTGGCC_L001' 'SECH_21-sech_PNF3_GTTTCG_L001')

i=$((${SLURM_ARRAY_TASK_ID} - 1))
python3 assemble_sech_genome.py ${sech_sample_ids[${i}]}

