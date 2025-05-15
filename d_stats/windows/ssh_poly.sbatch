#!/bin/bash
#SBATCH -p general
#sbatch -J ssh_poly
#SBATCH --cpus-per-task=1
#SBATCH -t 96:00:00
#SBATCH --mem=8g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/logs/ssh_poly_%A.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/logs/ssh_poly_%A.err


ml python/3.12.2 

###get polymorphic stats for windows across the genome for multiple populations
##python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_poly_windows.py 50000 mel \
##                                                                                          /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
##                                                                                          /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db


#get divergence stats for windows across the genome
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_poly_diff_windows.py 50000 mel \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

