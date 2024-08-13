#!/bin/bash
#SBATCH -p general
#sbatch -J hal2maf
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=4g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_maf_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_maf_%j.err


source /proj/matutelb/software/cactus-bin-v2.8.0/venv-cactus-v2.8.0/bin/activate
cactus-hal2maf --refGenome D_SIMULANS --targetGenomes D_SIMULANS,D_MELANOGASTER --dupeMode single --filterGapCausingDupes --chunkSize 500000 --batchCount 1 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/cactus_jobstore /proj/matutelb/data_share/drosophila_genomes/drosophila.hal /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_mel.maf
