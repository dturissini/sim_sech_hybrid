#!/bin/bash
#SBATCH -p general
#sbatch -J run_minimap2
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH --mem=8g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/inversion_check_minimap_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/assembly/logs/inversion_check_minimap_%j.err



module add minimap2/2.26
python3 inversion_check_minimap.py ${1} ${2}
