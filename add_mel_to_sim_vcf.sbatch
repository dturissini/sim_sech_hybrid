#!/bin/bash
#SBATCH -p general
#sbatch -J add_out_to_vcf
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00
#SBATCH --mem=4g
#SBATCH -o /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/add_outgroup_to_vcf_%j.out
#SBATCH -e /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/add_outgroup_to_vcf_%j.err


python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/add_outgroup_to_vcf.py D_SIMULANS D_MELANOGASTER /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/sim_mel.maf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/mel_aligned_sim.fasta /proj/matutelb/data_share/simulans_OA_resistance/simulans_sechellia.vcf.gz /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel.vcf.gz
