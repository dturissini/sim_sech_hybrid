chroms=("2L" "2R" "3L" "3R" "X")

for chrom in ${chroms[@]}; do
  job_file="/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_${chrom}/run_ripta_${chrom}.sh"
  
  echo "python3  /proj/matutelb/software/RIPTA/site_patterns.py -v /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_${chrom}.recode.vcf.gz \\" > ${job_file}
  echo "                                                        -m /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_sech_hybrid_ripta_meta.txt \\" >> ${job_file}
  echo "                                                        -f True \\" >> ${job_file}
  echo "                                                        -p1 simulans \\" >> ${job_file}
  echo "                                                        -p2 sim_sech_hybrid \\" >> ${job_file}
  echo "                                                        -p3 sechellia \\" >> ${job_file}
  echo "                                                        -p4 melanogaster \\" >> ${job_file}
  echo "                                                        -p /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_${chrom}/" >> ${job_file}

  echo "" >> ${job_file}
  
  echo "python3 /proj/matutelb/software/RIPTA/bootstrapping_one_thread.py -v /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel_${chrom}.recode.vcf.gz \\" >> ${job_file}
  echo "                                                                  -m /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/sim_sech_hybrid_ripta_meta.txt \\" >> ${job_file}
  echo "                                                                  -f \\" >> ${job_file}
  echo "                                                                  -p1 simulans \\" >> ${job_file}
  echo "                                                                  -p2 sim_sech_hybrid \\" >> ${job_file}
  echo "                                                                  -p3 sechellia \\" >> ${job_file}
  echo "                                                                  -p4 melanogaster \\" >> ${job_file}
  echo "                                                                  -cl 120905797 \\" >> ${job_file}
  echo "                                                                  -bs 100000 \\" >> ${job_file}
  echo "                                                                  -r 100 \\" >> ${job_file}
  echo "                                                                  -p /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/ripta_${chrom}/" >> ${job_file}

  chmod +x ${job_file}
done

