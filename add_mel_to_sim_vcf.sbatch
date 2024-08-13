#!/bin/bash
#SBATCH -p general
#sbatch -J add_out_to_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH --mem=4g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/add_outgroup_to_vcf_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/add_outgroup_to_vcf_%j.err


python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/add_outgroup_to_vcf.py D_SIMULANS D_MELANOGASTER /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_mel.maf /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/mel_aligned_sim.fasta /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.vcf.gz /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel.vcf.gz
