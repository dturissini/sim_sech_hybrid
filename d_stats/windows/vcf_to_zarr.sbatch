#!/bin/bash
#SBATCH -p general
#sbatch -J vcf_to_zarr
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=16g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/vcf_to_zarr_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/vcf_to_zarr_%j.err

python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/filtered_vcf_to_zarr.py 2L /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2L.recode.vcf.gz /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2L.zarr
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/filtered_vcf_to_zarr.py 2R /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2R.recode.vcf.gz /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_2R.zarr
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/filtered_vcf_to_zarr.py 3L /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3L.recode.vcf.gz /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3L.zarr
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/filtered_vcf_to_zarr.py 3R /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3R.recode.vcf.gz /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_3R.zarr
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/filtered_vcf_to_zarr.py X  /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_X.recode.vcf.gz  /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_X.zarr
