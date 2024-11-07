# sim_sech_hybrid background
Scripts for investigating introgression between Drosophila sechellia (sech) and D. simulans (sim) in the Seychelles islands. The analysis utilizes samples: D. sechellia from the Seychelles, cosmopolitan D. simulans from around the world, and sim_sech_hybrid (ssh) samples from the Seychelles. The sim_sech_hybrid flies are morphologically similar to D. simulans but show both male and female mating preference for sim_sech_hybrids against both D. simulans and D.sechellia flies. Sim_sech_hybrid flies also show low levels of resistance to octanoic acid which is found in Morinda citrifolia fruits, the primary food source for D. sechellia. Whereas, cosmopolitan D. simulans flies have no octanoic acid resistance. Sim_sech_hybrid flies, are therefore, thought to be D. simulans flies that have received introgression from D. sechellia in the Seycheeles. 

Previous attempts to identify introgressed tracts in sim_sech_hybrids have been unsuccessful (unpublished). It, therefore, appears that there are no large introgressed regions in the sim_sech_hybrid flies. This git repo contains scripts to better understand patterns of introgression focusing on genomic windows using the D+ statistic which performs better with small windows than the original ABA-BABA D statistic. Additionally, the analysis looks at levels of polymorphism in genomic windows in D. sechellia, a species known for very low levels of genomic polymorphism.

Commands run for most of the introgression steps described below can be found in sech_into_ssh_introgression_windows_steps.txt and sim_into_sech_introgression_windows_steps.txt. Additional text files for running the steps for other parts of the analysis can be found below.

Most of the scripts are written in python, R, and bash. Two sqlite databases are also created to store metadata and results simplifying processing and analyzing data. The analysis also makes extesnive use of submitting jobs to a slurm cluster.

## Data preparation
### Prepare vcf and zarr files
The starting vcf file only had samples from the D. simulans clade and lacked an outgroup. Alleles from the D. melanogaster reference genome were added to the vcf to serve as an outgroup. The vcf was then split by chromosomes and the chromomosome-specific vcfs then converted to the zarr file format to facilitate more efficiently obtaining genotype data for subsequenct steps. Commands for the steps to make the vcf and zarr files can be found in make_ssh_vcf_and_zarr_steps.txt.

1. Make two species alignment of D.simulans and D. melanogaster (mel) to facilitate identifying mel alleles based on sim genomic positions. 
    * cactus_hal2maf.sbatch usess the cactushal2maf program to read a  multispecies hal file and ouyt D. melabnogaster sequences aligned to the D. simulans reference genome.

2. Convert maf file to fasta file
    * maf_to_fasta.sbatch converts the maf file to a full genome (D. simulans) fasta file
    * Calls the python script make_outgroup_fasta_from_maf.py

3. Add D. melanogaster (mel) alleles to D. simulans clade vcf as an outgroup
    * add_mel_to_sim_vcf.sbatch processed the D. simulans clade vcf file and outputs a new vcf with the D. melanogaster alleles added. Sincewe're using the haploid D. melnogaster reference genome, all sites are assumed to be homozygous.
    * Calls the python script add_outgroup_to_vcf.py
    
4. Split vcf by chromosomal arm into separate vcf files
    * split_vcf_by_chr.sbatch uses vcftools to split the vcf file into spearate vcf files for each chromosome arm to allow for parallelization across chromosome arms by downstream analyses.

    
5. Convert vcf files to zarr files to speed up accessing genotype data by subsequent scripts
    * sbatch vcf_to_zarr.sbatch 
    * Calls the python script filtered_vcf_to_zarr.py


### Prepare metadata
The analysis relies on a number of sqlite database tables to store metadata related to the samples and populations. Steps to prepare the metadata can be found in make_ssh_metadata_steps.txt

1. Make chrom_lens table to store the length of chromosome arms for use when creating windows along chromosomes and populate it with a sql query

2. Make sample_species table to identify which species each sample is from. The table also contains the collection location for each sample (where available) and the relative column in the vcf file for that sample to simplify indexing when processing vcf and zarr files.
    * get_filtered_metadata.py is run to create and populate the table

3. Make sample_pop_link table to link each sample to populations (pops). Pops are either the the species (mel, sim, sech, ssh) or subpopulations defined by islands in the Seychelles. 
    * make_sample_pop_link.py is run to create and populate the table

4. Make pop_cols table assigning a specific color to each pop found in sample_pop_link to standardize pop colors across all plots and populate it with a sql query



## PCA
Commands for running the PCAs are found in pca/sim_sech_hybrid_mel_outgroup_pca_steps.txt. Briefly, a sqlite query is run on the linux command line to generate a text file with a subset of sample ids which is then used by plink to run the PCA. PCA plots are generated with ssh_mel_out_pca.R 


## D. sechellia into D. sechellia introgression
The sim_sech_hybrid (ssh) population of samples is suspected to be D. simulans flies that have received introgression from D. sechellia in the past. Collected in the Seychelles by Daniel Matute, the ssh flies show some degree of resistance to octanoic acid which is only seen in D. sechellia flies in the D. melanogaster clade. Given sympatry between the populations, introgression seems likely. The steps described below attempt to identify and quanitfy historic introgression from D. sechellia into ssh flies.

### D stats - genomic
D and D+ statistics were calculated for each chromosome arm using RIPTA. 

1. make_ripta_job_files.sh creates chromosome arm specific bash scripts that can be subsequently submitted as slurm jobs. Each bash script has two steps: 1) run the RIPTA site_patterns.py program to generate counts of the ABBA BABA site patterns and 2) the RIPTA bootstrapping_one_thread.py for performing bootstrapping to get distributions of ABBA BABA site counts.

2. run_ripta.sbatch runs the bash scripts created in teh previous step

3. ssh_dplus_bootstrap.R plots the D and D+ values as well as the bootstrap distributions

### Genomic windows
After identifying evidence of introgression at the genomic level, we next applied a window-based approach to try to identify specific regions with introgression. The D+ was statistic was used since it's been shown to have much better accuracy than D when applioed to genomic windows. Commands for running the steps described below can be found in sech_into_ssh_introgression_windows_steps.txt

Note that this section of the analysis consistently used sim, ssh, and mel as populations 1, 2, and 4 respectively within the ABBA-BABA D statistic framework, but multiple populations of D. sechellia were tested as population 3 including all sech samples as well as subpopulations defined based on geography of the Seychelles archipelago. Samples for the subpopulations were recorded in the sample_pop_link database table described in the Prepare Metadata section above.

#### D stats
1. get_d_stat_windows.py reads zarr files and divides the genome into windows and records D statistics and ABBA, BABA, ABAA, and BAAA counts in a database table with the prefix "d_stat_win_". 

2. d_plus_windows.R reads the D stat and ABBA BABA count data from the database table and makes a series of plots in a pdf file.
    * Distribution of number of polymorphic sites per window to help set useful cutoffs to remove windows with too few sites
    * Histograms of ABBA BABA counts to help visualize any outliers and set a useful threshold for windows with too informative sites
    * Scatterplot of nuymber of sites versus total AB sites with red lines showing cutoff thresholds
    * D+ distribution with a red line denoting the 99th quantile used for identifying outliers in downstream analyses
    * Scatterplots of different AB counts versus D+
    * Traces of per-window D+ across each chromosome arm with significant SNPs from a GWAS run on octanoic acid resistance in ssh flies
    * Scatterplots of D+ versus the component ABBA-BABA and BAAA-ABAA components to help better understand the statistic and visualize outliers
    * Scatterplot of GWAS p-values vs D+ for the respective genomic window for each significant SNP
    * Scatterplots of D+ versus the percent of sites for the window that are ABBA or BAAA to help better understand their respective contributions 


#### Polymorphism
Polymorphism was then measured both within and between populations for genomic windows and compared against the D+ results. D. sechellia, an island endemic species, has famously low polymorphism.

1. get_poly_windows.py reads zarr files and calculates pi and the derived allele frequency for each pop present in the sample_pop_link database table. Results are stored in a database table with the prefix "poly_win_". Additional, the site frequency spectrum (SFS) is calculated for each pop with results stroed in a table with the prefix "sfs_der_freq_".

2. get_poly_diff_windows.py reads zarr files and calculates several between population measurements of variation (Dxy, Fst, and difference of derived allele frequencies). Between population comparisons are only run on the most relevant subset of combinations of pops to reduce unnecessary computation. Results are stored in a database table with the prefix "poly_diff_win_".

3. poly_windows.R generates a pdf with plots for both polymorphism and D stat results. The R script takes as input the window size but also the four populations used for D+ statistic. The 3 pops (excluding outgroup) whose polymorphism data is plotted are obtained by parsing this string (e.g. sim_ssh_sech_mel).
    * Distribution of number of polymorphic sites per window to help set useful cutoffs to remove windows with too few sites
    * Histograms of per-window pi for sim, sech, and ssh
    * Histograms of per-window derived allele frequencies (mel as outgroup) for sim, sech, and ssh
    * Histograms of per-window pairwise differences for derived allele frequencies (mel as outgroup) for sim, sech, and ssh
    * Histograms of Dxy for all pariwise combinations of sim, sech, and ssh
    * Histograms of Fst for all pairwise combinations of sim, sech, and ssh
    * Per site derived allele frequency distributions for sim, sech, and ssh
    * Traces of per-window values for several measurements to help visualize both their distribtuons across the genome and any potential correlations.
        * Significant GWAS SNPs and their p-values for octanoic acid resistance in ssh
        * AB site counts
        * D+
        * pi
        * Dxy
        * Fst
        * Derived allele frequencies
    * If sech is one of the pops, then a number of additional plots are made
        * Scatterplot of pi sech versus D+ is plotted to see if D+ outlier windows have higher pi sech
        * Scatterplot of BAAA-ABAA vs pi sech with points colored to reflect D+ value
        * Scatterplot of BAAA-ABAA vs D+ with points colored and size-scaled to reflect pi sech value
        * Scatterplot of ABBA-BABA vs D+ with points colored and size-scaled to reflect pi sech value
      

#### Outlier windows
Outlier windows for pi sech, D+  were identified using the 99th quantile as the cutoff, and the total number of alleles and number of derived alleles for each pop was recorded for each polymorphic site within these windows. Results were stored in a database table. The outgroup and derived alleles were also recorded for each site in another database table. Additionally 50 random windows were selected and the same values recorded to serve as null distributions. 

  1. get_outlier_win_sites.py takes an input parameter identifying whether outlier pi sech, outlier D+, or random windows should be identified and processed.

  2. get_outlier_win_neighbor_sites records the same information as get_outlier_win_sites.py but processes windows neighboring outlier windows. This script is run for pi sech outlier windows to record information for snps which may still be part of the diverged haplotypes (observed in subsequent steps) that could extend beyon the boundaries of the outlier windows.

  3. outlier_wins.R generates a pdf with a page for each of the outlier windows: pi sech, D+, or random which are passed in as a parameter. Each window plot has several panels:
      * Scatterplots of dervied allele frequencies for each polymorphic site in the window. There is a separate panel for each of sim, sech, and ssh
      * A histogram to the right of the scatterplot showing the distribution of derived allele frequencies.
      * A separate scatterplot of derived allele frequencies for all of the sech pops (defined by Seychelles geography)
      * Traces of histogram counts showing the distribution of derived allele frequencies for each sech pop.
      * Additionally, transparent gray bars show the location of annotated genes in the window and vertical black lines indicate significant SNPs from the ssh octanoic acid GWAS. D+ and the genomic D+ quantile are also printed in the upper-right corner.

  4. get_outlier_win_alleles.py records the number of derived alleles for each site for each sample for all polymorphic sites in out lier windows (pi sech, D+, or random) and stores the results in a database table
  
  5. make_outlier_pi_sech_win_adj_allele_dists.py seeks to see how much each sample differs from the sech anro pop in outlier windows (pi sech, D+, or random) by summing the number of alleles that differ from sech anro as a rough distance metric. It takes advantage of the extremely low polymorphism in D. sechellia, and as a bit of a shortcut compares against a single sech anro sample as a reference. sech anro was chosen as the reference since it is low polymorphism and has consistent haplotypes in outlier pi sech windows that often differ from sech denis and sech praslin pops. The distance metric is helpful for easily visualizing the divergent haplotypes seen in sech for many of the outlier D+ and pi sech windows.
  
  6. outlier_sech_alleles.R makes a pdf with a page for each outlier window (pi sech, D+, or random). The plots contain two panels. The top panel contains the same pi sech derived allele frequency scatterplot as in outlier_wins.R. The bottom panel shows per site alleles for each sech sample on its own line. The alleles are colored by their similarity to a reference sech anro sample with blue being homozygous anro, green bing heterozygous, and red being homozygous non-sech anro. The plots are used to visualize the unexpectedly differentiated haplotypes we found in sech samples. Colored boxes on the left of the plot show how different each sample is from sech anro for the outlier pi sech alleles using the metric calculated in make_outlier_pi_sech_win_adj_allele_dists.py as described above.
  
  7. get_sech_anro_dxy_windows.py calcualtes Dxy between each sech pop and sech anro for windows across the genome and stores the results in a database table.
  
  8. sech_anro_dxy_windows.R generates a pdf with plots showing the window sech anro Dxy results from the previous step. 
      * Scatterplot of sech anro Dxy (averaged over all windows) against the sech anro allele distance calculated in step 5 above. 
      * Scatterplots similar to above but removing outlier pi sech windows above different thresholds to see if genomic patterns still persist or are driven by outlier windows.
      * Traces of window Dxy values across each chromosome arm for each sech pop compared against sech anro
      * Histograms for each sech sample of Dxy compared against sech anro
  
  
## D. simulans into D. sechellia introgression
The initial focus of the analysis was to look for introgression from D. sechellia into D. simulans resulting in the ssh population of flies. However, a very negative peak for D+ in the middle of chromosome 3R suggested there had also been gene flow from D. simulans into D. sechellia. sim_into_sech_introgression_windows_steps.txt contains the steps for this analysis. Here we split D. sechellia into two populations based on the population structure oberved in the D. sechellia into D. simulans analysis. sech depr contains all samples from Denis and Praslin and sech base the remaining samples (except those from an unknown location). sech base was used as the recipient population for introgression from D. simulans based on observed results. sech depr, sech denis, and sech pras were all run as the other D.sechellia population. 

Commands for running the steps described below can be found in sim_into_sech_introgression_windows_steps.txt are mostly the same as those described above for D. sechellia into D. Sinulans introgression with one additional R script:

  1. poly_windows_comp.R creates a pdf with one page for each chromosome arm comparing ABBA BABA counts and D+ for sim_ssh_sech_mel and sechdepr_sech_base_sim_mel.An additional panel showed avergae pi sech for each window for sech, (sech depr, sech denis, or sech pras), and sech base. Transparent grey bars show pi sech outlier windows (top 1%). These plots are helpful for seeing that outlier pi sech windows overlap with outlier D+ windows for introgression between D. sechellia and D. simulans in both directions. The R script also runs hypergeometric tests to show a significnat overlap between pi sech outliers and D+



## Inversion checks
To investigate putative inversions, several D. sechellia lines were sequenced using Oxford Nanopore, and their genomes were assembled using Canu. Note that the fly lines were intially sequenced over a decade ago, anbd many were, unfortunately, lost during COVID. So it was not always possible to sequence the specif D. sechellia lines from Denis, and additional Denis lines were also sequenced with the expectation that any inversions segregating in the population could also be found in them too.

Note that the inversion checks are actively being worked on, and we are currently waiting on the results of the next round of Oxford Nanopore sequencing following the lab updating sample preparation protocols.

### Genome assembly
CLEANUP sech assembly steps txt

### Look for putative inversions
inversion_check_steps.txt

#### Raw reads

##### Map long reads to reference genome
python3 process_inversion_overlaps.py trimmedReads

sbatch inversion_check_lastz_reads.sbatch sech_Denis_NF69

identify sech Anro alleles for polymorphic sites within putative inversions by looking at illumina data
python3 get_inversion_sech_alleles.py /work/users/d/t/dturissi/drosophila/ssh/introgression/sim_sech_outgroup_mel \

get number of derived and sech anro alleles for reads overlapping putative inversions
python3 get_derived_alleles_from_mafs_reads.py 

look at anro allele concordance for reads overlapping putative inversions
Rscript inversion_check_per_der_reads.R


look at anro allele concordance for reads overlapping putative inversions, coloring reads by their amount of anro allele concordance
Rscript inversion_check_per_der_reads_inversion_plots.R




##### Plot aligned reads that overlap putative inversions
Rscript inversion_check_contig_map.R 50000 sech_Denis_NF22 trimmedReads


#### Assembled contigs

##### Map contigs to reference genome
python3 process_inversion_overlaps.py contigs

LastZ used to get sequnece alignments
sbatch inversion_check_lastz.sbatch sech_Denis_NF22

get number of derived and sech anro alleles  for contigs overlapping putative inversions
python3 get_derived_alleles_from_contig_mafs.py 


##### Plot aligned contigs that overlap putative inversions
Rscript inversion_check_contig_map.R 50000 sech_Denis_NF22 contigs


