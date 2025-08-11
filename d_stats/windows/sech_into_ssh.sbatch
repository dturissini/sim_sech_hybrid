#!/bin/bash
#SBATCH -p general
#sbatch -J sech_into_ssh
#SBATCH --cpus-per-task=1
#SBATCH -t 240:00:00
#SBATCH --mem=16g
#SBATCH -o /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/logs/sech_into_ssh_%j.out
#SBATCH -e /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/logs/sech_into_ssh_%j.err


ml python/3.12.2 

                                                                                                       
#get window D statistics for sim_ssh_sech_mel using all sech samples
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_d_stat_windows.py sim \
                                                                                            ssh \
                                                                                            sech \
                                                                                            mel \
                                                                                            50000 \
                                                                                            /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                            /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db



#get window D statistics for sim_ssh_sechbase_mel using only sech samples from Anro, La Digue, and Marianne
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_d_stat_windows.py sim \
                                                                                            ssh \
                                                                                            sechbase \
                                                                                            mel \
                                                                                            50000 \
                                                                                            /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                            /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db
                                                                                                       


#get window D statistics for sim_ssh_sechpras_mel using only sech samples from Praslin
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_d_stat_windows.py sim \
                                                                                            ssh \
                                                                                            sechpras \
                                                                                            mel \
                                                                                            50000 \
                                                                                            /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                            /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db
                                                                                            
                                                                                                                                                                                                              
#create pdf reports for D statistic windows                                                                                                       
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/d_plus_windows.R 50000 sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/d_plus_windows.R 50000 sim_ssh_sechbase_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/d_plus_windows.R 50000 sim_ssh_sechpras_mel




#batch file to run mulktiple steps
##sbatch /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_poly.sbatch


#get polymorphic stats for windows across the genome for multiple populations
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_poly_windows.py 50000 mel \
                                                                                          /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                          /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db


#get divergence stats for windows across the genome
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_poly_diff_windows.py 50000 mel \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db



#create pdfs containing plots showing polymorphism, D statistic, and divergence stats across the genome
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/poly_windows.R 50000 sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/poly_windows.R 50000 sim_ssh_sechbase_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/poly_windows.R 50000 sim_ssh_sechpras_mel



#Identify outlier windows with high or low D+ values and record the number of derived and total alleles for each pop for all sites within these windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_sites.py 50000 \
                                                                                               d_plus \
                                                                                               sim_ssh_sech_mel \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

#Randomly select 50 windows to serve as an empirical null distribution and record the number of derived and total alleles for each pop for all sites within these windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_sites.py 50000 \
                                                                                               random \
                                                                                               sim_ssh_sech_mel \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

#Identify outlier windows with high pi in sechellia and record the number of derived and total alleles for each pop for all sites within these windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_sites.py 50000 \
                                                                                               pi_sech \
                                                                                               sim_ssh_sech_mel \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                               /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db


#make plots for each outlier window showing the derived allele frequencies for simulans, sechellia, and sim_sech_hybrids
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_wins.R 50000 d_plus sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_wins.R 50000 random sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_wins.R 50000 pi_sech sim_ssh_sech_mel


#additional steps to get bordering windows for pi_sech peaks to help better identify breakpoints for putative inversions
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_neighbor_sites.py 50000 \
                                                                                                        pi_sech \
                                                                                                        sim_ssh_sech_mel \
                                                                                                        /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                        /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_wins.R 50000 neighbor_pi_sech sim_ssh_sech_mel



#get alleles for all sech samples for sites in outlier D+ windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_alleles.py 50000 \
                                                                                                 sech \
                                                                                                 d_plus \
                                                                                                 sim_ssh_sech_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

################
#debug


#prep
sqlite3 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/debug.db
attach database '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db' as s;

create table sample_species as select * from s.sample_species;
create index idx_ss_sid on sample_species(sample_id);

create table sample_pop_link as select * from s.sample_pop_link;
create index idx_spl_pop on sample_pop_link(pop);
create index idx_spl_sample_id on sample_pop_link(sample_id);

create table d_stat_win_50000_sim_ssh_sech_mel as select * from s.d_stat_win_50000_sim_ssh_sech_mel;
create index idx_d_swe on d_stat_win_50000_sim_ssh_sech_mel(start,end);
create index idx_d_see on d_stat_win_50000_sim_ssh_sech_mel(end);


create table outlier_d_plus_win_sites_50000_sim_ssh_sech_mel as select * from s.outlier_d_plus_win_sites_50000_sim_ssh_sech_mel;

create index idx_odwps5_pos_50000d_plus on outlier_d_plus_win_sites_50000_sim_ssh_sech_mel(pos);
create index idx_odwps5_win_id_50000d_plus on outlier_d_plus_win_sites_50000_sim_ssh_sech_mel(win_id);
create index idx_odwps5_pop_50000d_plus on outlier_d_plus_win_sites_50000_sim_ssh_sech_mel(pop);


python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_alleles.py 50000 \
                                                                                                 sech \
                                                                                                 d_plus \
                                                                                                 sim_ssh_sech_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/debug.db


sqlite3 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/debug.db
attach database 'ssh_d_win_old.db' as o;

select count(*), sum(case when n.num_der_alleles != o.num_der_alleles then 1 else 0 end)
from outlier_d_plus_win_alleles_sech_50000_sim_ssh_sech_mel n, o.outlier_d_plus_win_alleles_sech_50000_sim_ssh_sech_mel o
where n.chrom = o.chrom
and n.pos = o.pos
and n.sample_id = o.sample_id;
446292|0



sbatch debug_vcf_allele_diffs.sbatch
sacct -j 10044527

Rscript debug_vcf_avg_depths.R

#################      

#get alleles for all sech samples for sites in the random "outlier" windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_alleles.py 50000 \
                                                                                                 sech \
                                                                                                 random \
                                                                                                 sim_ssh_sech_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db


#get alleles for all sech samples for sites in outlier pi sech windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_alleles.py 50000 \
                                                                                                 sech \
                                                                                                 pi_sech \
                                                                                                 sim_ssh_sech_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db


#look at sech and ssh alleles to see if putative sech inversions are also present in ssh (spoiler: they don't seem to be)
#get alleles for all ssh samples for sites in outlier D+ windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_alleles.py 50000 \
                                                                                                 ssh \
                                                                                                 d_plus \
                                                                                                 sim_ssh_sech_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db


#get alleles for all ssh samples for sites in outlier pi sech windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_outlier_win_alleles.py 50000 \
                                                                                                 ssh \
                                                                                                 pi_sech \
                                                                                                 sim_ssh_sech_mel \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                 /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db



#Get the allele distance (sum of numnber of differences) between each sample and a sech anro sample (SECH_3-sech_Anro_B3_TTAGGC_L001). 
#This measure takes advantage of the very low polymorphism in sechellia where there are very few polymorphic sites in sech anro.
#The distance is helpful for easily visualizing the divergent haplotypes seen in sech for many of the outlier D+ and pi sech windows.
#get allele distances for D+ outlier windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/make_outlier_pi_sech_win_adj_allele_dists.py 50000 \
                                                                                                                   d_plus \
                                                                                                                   sim_ssh_sech_mel \
                                                                                                                   /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

#get allele distances for random "outlier" windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/make_outlier_pi_sech_win_adj_allele_dists.py 50000 \
                                                                                                                   random \
                                                                                                                   sim_ssh_sech_mel \
                                                                                                                   /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

#get allele distances for pi sech outlier windows
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/make_outlier_pi_sech_win_adj_allele_dists.py 50000 \
                                                                                                                   pi_sech \
                                                                                                                   sim_ssh_sech_mel \
                                                                                                                   /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db





#plot the alleles for each sample for outlier windows.
#Alleles are colored by whether they have a sech anro allele or not using SECH_3-sech_Anro_B3_TTAGGC_L001 as the reference for polarizing alleles.
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_sech_alleles.R 50000 d_plus sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_sech_alleles.R 50000 random sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_sech_alleles.R 50000 pi_sech sim_ssh_sech_mel






Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_pi_sech_win_adj_allele_dists.R 50000 sim_ssh_sech_mel
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/outlier_sech_ssh_alleles.R 50000 sim_ssh_sech_mel


#Calculate dxy between samples and sech anro for windows
#These dxy values can be used to see if samples that are more or less diverged to sech anro in outlier windows are also similarly diverged genome wide
python3 /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/get_sech_anro_dxy_windows.py 50000 \
                                                                                                   /work/users/d/t/dturissi/drosophila/ssh/genotyping/zarrs/sim_sech_hybrid_filtered_biallelic_snp_mel_outgroup \
                                                                                                   /work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db

#Make pdf showing sech anro dxy values for sech and ssh pops across the genome and how they relate to pi sech
Rscript /work/users/d/t/dturissi/drosophila/ssh/introgression/scripts/sech_anro_dxy_windows.R 50000 sim_ssh_sech_mel


