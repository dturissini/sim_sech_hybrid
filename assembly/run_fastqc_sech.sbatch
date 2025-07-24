#!/bin/bash
#SBATCH -p general
#sbatch -J fastqc_sech
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH --mem=8g
#SBATCH --array=1-9
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/run_fastqc_sech_%A_%a.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/run_fastqc_sech_%A_%a.err


sech_sample_ids=('sech_Denis_NF22' 'sech_Denis_NF69' 'sech_Denis_NF712' 'sech_Denis_papaya152' 'sech_Anro_114_female' 'sech_Denis_NF45_female' 'DenisNF134' 'SECH_20-sech_PNF5_GTGGCC_L001' 'SECH_21-sech_PNF3_GTTTCG_L001')


ml fastqc 

i=$((${SLURM_ARRAY_TASK_ID} - 1))

fastqc -f fastq --memory 7500 -o fastqc /work/users/d/t/dturissi/drosophila/ssh/assembly/fastqs/${sech_sample_ids[${i}]}.fastq

