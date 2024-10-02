#!/bin/bash
#SBATCH -p general
#sbatch -J run_lastz
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH --mem=16g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/inversion_check_lastz_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/inversion_check_lastz_%j.err



module add lastz/1.03.73   
python3 inversion_check_lastz.py ${1}
