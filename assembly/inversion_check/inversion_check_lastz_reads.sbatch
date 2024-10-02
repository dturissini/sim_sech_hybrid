#!/bin/bash
#SBATCH -p general
#sbatch -J run_lastz
#SBATCH --cpus-per-task=1
#SBATCH -t 96:00:00
#SBATCH --mem=16g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/inversion_check_lastz_reads_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/inversion_check_lastz_reads_%j.err



module add lastz/1.03.73   
python3 inversion_check_lastz_reads.py ${1}

python3 get_derived_alleles_from_read_mafs.py