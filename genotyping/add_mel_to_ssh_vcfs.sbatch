#!/bin/bash
#SBATCH -p general
#sbatch -J add_mel_to_ssh_vcf
#SBATCH -t 240:00:00
#SBATCH --mem=4g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/genotyping/logs/add_mel_to_ssh_vcf_%A_%a.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/genotyping/logs/add_mel_to_ssh_vcf_%A_%a.err


python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/add_outgroup_to_vcf.py D_SIMULANS D_MELANOGASTER \
                                                                                             /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/mel_aligned_sim.fasta \
                                                                                             /work/users/d/t/dturissi/drosophila/ssh/genotyping/vcfs/sim_sech_hybrid.vcf.gz \
                                                                                             /work/users/d/t/dturissi/drosophila/ssh/genotyping/vcfs/sim_sech_hybrid_mel_outgroup.vcf.gz
