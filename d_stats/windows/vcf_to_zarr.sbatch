#!/bin/bash
#SBATCH -p general
#sbatch -J vcf_to_zarr
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=16g
#SBATCH -o /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/vcf_to_zarr_%j.out
#SBATCH -e /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/vcf_to_zarr_%j.err

python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/filtered_vcf_to_zarr.py 2L /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2L.recode.vcf.gz /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2L.zarr
python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/filtered_vcf_to_zarr.py 2R /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2R.recode.vcf.gz /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_2R.zarr
python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/filtered_vcf_to_zarr.py 3L /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3L.recode.vcf.gz /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3L.zarr
python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/filtered_vcf_to_zarr.py 3R /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3R.recode.vcf.gz /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_3R.zarr
python3 /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/scripts/filtered_vcf_to_zarr.py X  /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_X.recode.vcf.gz  /proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/sim_sech_outgroup_mel_X.zarr
