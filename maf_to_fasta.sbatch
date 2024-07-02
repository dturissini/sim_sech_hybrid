#!/bin/bash
#SBATCH -p general
#sbatch -J maf2fasta
#SBATCH --cpus-per-task=1
#SBATCH -t 4:00:00
#SBATCH --mem=4g
#SBATCH -o /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/maf_to_fasta_%j.out
#SBATCH -e /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/maf_to_fasta_%j.err


python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/make_outgroup_fasta_from_maf.py D_SIMULANS D_MELANOGASTER /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/sim_mel.maf /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/mel_aligned_sim.fasta 
