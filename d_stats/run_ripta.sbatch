#!/bin/bash
#SBATCH -p general
#sbatch -J ripta
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=16g
#SBATCH --array=1-5
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/run_ripta_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/run_ripta_%j.err

/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_2L/run_ripta_2L.sh
/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_2R/run_ripta_2R.sh
/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_3L/run_ripta_3L.sh
/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_3R/run_ripta_3R.sh
/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_X/run_ripta_X.sh