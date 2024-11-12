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
        * Inputs: 
            * reference species name
            * outgroup species name
            * maf file path
            * reference genome fasta file path

3. Add D. melanogaster (mel) alleles to D. simulans clade vcf as an outgroup
    * add_mel_to_sim_vcf.sbatch processed the D. simulans clade vcf file and outputs a new vcf with the D. melanogaster alleles added. Sincewe're using the haploid D. melnogaster reference genome, all sites are assumed to be homozygous.
    * Calls the python script add_outgroup_to_vcf.py
        * Inputs: 
            * reference species name
            * outgroup species name
            * maf file path
            * outgroup genome fasta file path
            * input vcf file path
            * output vcf file path
    
4. Split vcf by chromosomal arm into separate vcf files
    * split_vcf_by_chr.sbatch uses vcftools to split the vcf file into spearate vcf files for each chromosome arm to allow for parallelization across chromosome arms by downstream analyses.

    
5. Convert vcf files to zarr files to speed up accessing genotype data by subsequent scripts
    * sbatch vcf_to_zarr.sbatch 
    * Calls the python script filtered_vcf_to_zarr.py
        * Inputs: 
            * chromosome name
            * vcf file path
            * output zarr file path


### Prepare metadata
The analysis relies on a number of sqlite database tables to store metadata related to the samples and populations. ssh_d_win.db is created to store the metadata. Steps to prepare the metadata can be found in make_ssh_metadata_steps.txt

1. Make chrom_lens table to store the length of chromosome arms for use when creating windows along chromosomes and populate it with a sql query

2. Make sample_species table to identify which species each sample is from. The table also contains the collection location for each sample (where available) and the relative column in the vcf file for that sample to simplify indexing when processing vcf and zarr files.
    * get_filtered_metadata.py is run to create and populate the table
        * Inputs: 
            * tab-delimited metadata file path
            * vcf file path
            * sqlite database file path (note sqlite will make the file if it doesn't already exist)

3. Make sample_pop_link table to link each sample to populations (pops). Pops are either the the species (mel, sim, sech, ssh) or subpopulations defined by islands in the Seychelles. 
    * make_sample_pop_link.py is run to create and populate the table
        * Inputs: 
            * tab-delimited pop file path with two columns: pop name, species name
            * sqlite database file path


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

The same ssh_d_win.db database that contains the metadata is also used to restore results for the analyses described below.

Note that this section of the analysis consistently used sim, ssh, and mel as populations 1, 2, and 4 respectively within the ABBA-BABA D statistic framework, but multiple populations of D. sechellia were tested as population 3 including all sech samples as well as subpopulations defined based on geography of the Seychelles archipelago. Samples for the subpopulations were recorded in the sample_pop_link database table described in the Prepare Metadata section above.

#### D stats
1. get_d_stat_windows.py reads zarr files and divides the genome into windows and records D statistics and ABBA, BABA, ABAA, and BAAA counts in a database table with the prefix "d_stat_win_". 
    * Inputs: 
        * D stat population 1 name
        * D stat population 2 name
        * D stat population 3 name
        * D stat population 4 name
        * window size
        * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
        * sqlite database file path

2. d_plus_windows.R reads the D stat and ABBA BABA count data from the database table and makes a series of plots in a pdf file.
    * Inputs: 
        * window size
        * D stat pop string (e.g. sim_ssh_sech_mel)
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
    * Inputs: 
        * window size
        * outgroup species pop
        * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
        * sqlite database file path

2. get_poly_diff_windows.py reads zarr files and calculates several between population measurements of variation (Dxy, Fst, and difference of derived allele frequencies). Between population comparisons are only run on the most relevant subset of combinations of pops to reduce unnecessary computation. Results are stored in a database table with the prefix "poly_diff_win_".
    * Inputs: 
        * window size
        * outgroup species pop
        * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
        * sqlite database file path

3. poly_windows.R generates a pdf with plots for both polymorphism and D stat results. The R script takes as input the window size but also the four populations used for D+ statistic. The 3 pops (excluding outgroup) whose polymorphism data is plotted are obtained by parsing this string (e.g. sim_ssh_sech_mel).
    * Inputs: 
        * window size
        * D stat pop string (e.g. sim_ssh_sech_mel)
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
      * Inputs: 
          * window size
          * outlier type (d_plus, pi_sech, random)
          * D stat pop string (e.g. sim_ssh_sech_mel)
          * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
          * sqlite database file path

  2. get_outlier_win_neighbor_sites.py records the same information as get_outlier_win_sites.py but processes windows neighboring outlier windows. This script is run for pi sech outlier windows to record information for snps which may still be part of the diverged haplotypes (observed in subsequent steps) that could extend beyon the boundaries of the outlier windows.
      * Inputs: 
          * window size
          * outlier type (d_plus, pi_sech, random)
          * D stat pop string (e.g. sim_ssh_sech_mel)
          * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
          * sqlite database file path

  3. outlier_wins.R generates a pdf with a page for each of the outlier windows: pi sech, D+, or random which are passed in as a parameter. Each window plot has several panels:
      * Inputs: 
          * window size
          * outlier type (d_plus, pi_sech, random)
          * D stat pop string (e.g. sim_ssh_sech_mel)
      * Scatterplots of dervied allele frequencies for each polymorphic site in the window. There is a separate panel for each of sim, sech, and ssh
      * A histogram to the right of the scatterplot showing the distribution of derived allele frequencies.
      * A separate scatterplot of derived allele frequencies for all of the sech pops (defined by Seychelles geography)
      * Traces of histogram counts showing the distribution of derived allele frequencies for each sech pop.
      * Additionally, transparent gray bars show the location of annotated genes in the window and vertical black lines indicate significant SNPs from the ssh octanoic acid GWAS. D+ and the genomic D+ quantile are also printed in the upper-right corner.

  4. get_outlier_win_alleles.py records the number of derived alleles for each site for each sample for all polymorphic sites in out lier windows (pi sech, D+, or random) and stores the results in a database table
      * Inputs: 
          * window size
          * pop to record alleles for (e.g. sech, ssh)
          * outlier type (d_plus, pi_sech, random)
          * D stat pop string (e.g. sim_ssh_sech_mel)
          * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
          * sqlite database file path
  
  5. make_outlier_pi_sech_win_adj_allele_dists.py seeks to see how much each sample differs from the sech anro pop in outlier windows (pi sech, D+, or random) by summing the number of alleles that differ from sech anro as a rough distance metric. It takes advantage of the extremely low polymorphism in D. sechellia, and as a bit of a shortcut compares against a single sech anro sample as a reference. sech anro was chosen as the reference since it is low polymorphism and has consistent haplotypes in outlier pi sech windows that often differ from sech denis and sech praslin pops. The distance metric is helpful for easily visualizing the divergent haplotypes seen in sech for many of the outlier D+ and pi sech windows.
      * Inputs: 
          * window size
          * outlier type (d_plus, pi_sech, random)
          * D stat pop string (e.g. sim_ssh_sech_mel)
          * sqlite database file path
  
  6. outlier_sech_alleles.R makes a pdf with a page for each outlier window (pi sech, D+, or random). The plots contain two panels. The top panel contains the same pi sech derived allele frequency scatterplot as in outlier_wins.R. The bottom panel shows per site alleles for each sech sample on its own line. The alleles are colored by their similarity to a reference sech anro sample with blue being homozygous anro, green bing heterozygous, and red being homozygous non-sech anro. The plots are used to visualize the unexpectedly differentiated haplotypes we found in sech samples. Colored boxes on the left of the plot show how different each sample is from sech anro for the outlier pi sech alleles using the metric calculated in make_outlier_pi_sech_win_adj_allele_dists.py as described above.
      * Inputs: 
          * window size
          * outlier type (d_plus, pi_sech, random)
          * D stat pop string (e.g. sim_ssh_sech_mel)
  
  7. get_sech_anro_dxy_windows.py calcualtes Dxy between each sech pop and sech anro for windows across the genome and stores the results in a database table.
      * Inputs: 
          * window size
          * prefix for zarr file path (file name minus chromosome and extension, e.g. full path minus '_2L.zarr')
          * sqlite database file path
  
  8. sech_anro_dxy_windows.R generates a pdf with plots showing the window sech anro Dxy results from the previous step. 
      * Inputs: 
          * window size
          * D stat pop string (e.g. sim_ssh_sech_mel)
      * Scatterplot of sech anro Dxy (averaged over all windows) against the sech anro allele distance calculated in step 5 above. 
      * Scatterplots similar to above but removing outlier pi sech windows above different thresholds to see if genomic patterns still persist or are driven by outlier windows.
      * Traces of window Dxy values across each chromosome arm for each sech pop compared against sech anro
      * Histograms for each sech sample of Dxy compared against sech anro
  
  
## D. simulans into D. sechellia introgression
The initial focus of the analysis was to look for introgression from D. sechellia into D. simulans resulting in the ssh population of flies. However, a very negative peak for D+ in the middle of chromosome 3R suggested there had also been gene flow from D. simulans into D. sechellia. sim_into_sech_introgression_windows_steps.txt contains the steps for this analysis. Here we split D. sechellia into two populations based on the population structure oberved in the D. sechellia into D. simulans analysis. sech depr contains all samples from Denis and Praslin and sech base the remaining samples (except those from an unknown location). sech base was used as the recipient population for introgression from D. simulans based on observed results. sech depr, sech denis, and sech pras were all run as the other D.sechellia population. 

Commands for running the steps described below can be found in sim_into_sech_introgression_windows_steps.txt are mostly the same as those described above for D. sechellia into D. Sinulans introgression with one additional R script:

  1. poly_windows_comp.R creates a pdf with one page for each chromosome arm comparing ABBA BABA counts and D+ for sim_ssh_sech_mel and sechdepr_sech_base_sim_mel.An additional panel showed avergae pi sech for each window for sech, (sech depr, sech denis, or sech pras), and sech base. Transparent grey bars show pi sech outlier windows (top 1%). These plots are helpful for seeing that outlier pi sech windows overlap with outlier D+ windows for introgression between D. sechellia and D. simulans in both directions. The R script also runs hypergeometric tests to show a significnat overlap between pi sech outliers and D+
      * Inputs: 
          * window_size
          * First D stat pop string (e.g. sim_ssh_sech_mel)
          * Second D stat pop string (e.g. sim_ssh_sech_mel)


## Inversion checks
To investigate putative inversions, several D. sechellia lines were sequenced using Oxford Nanopore, and their genomes were assembled using Canu. Note that the fly lines were intially sequenced over a decade ago, anbd many were, unfortunately, lost during COVID when access to the lab was restricted. So it was not always possible to sequence the specif D. sechellia lines from Denis, and additional Denis lines were also sequenced with the expectation that any inversions segregating in the population could also be found in them too.

Note that the inversion checks are actively being worked on, and we are currently waiting on the results of the next round of Oxford Nanopore sequencing following the lab updating sample preparation protocols.

### Genome assembly
Steps used to assemble the genomes using canu can be found in sech_assembly_steps.txt.

  1. Install canu
  2. Create directories for results
  3. Prepare fastq files
  4. Run canu for each sample toa ssemble the genome by submitting a slurm array using canu_sech_genomes.sbatch. Note that several of the canu runs expectedly stopped, but were able to restart by just resubmitting the can job for the affected samples.
  5. Download the D. simulans reference genome and rename chromosomes (This step could also be moved to the inversion check stage)

### Look for putative inversions
inversion_check_steps.txt contains the steps run to look for putative inversions identifed from the diverged haplotypes in sech denis and desc praslin seen in the introgression analysis. In many cases neighobroing windows were combined, and the plots were manually reviewed to define rough start and end positions where the breakpoints for possible inversions should be found. Inversion checks were run for both raw reads and assembled contigs. The focus was to try to identify a read or contig that was split and mapped to two different locations to the reference D. simulans genome indicating an inversion breakpoint within the read/contig. Since many of the original D. sechellia lines had died off in the lab and were no longer available, other sech denis samples were sequenced, and their alleles were compared to the reference sech anro sample to try to identify assembled genomes that had the diverged haplotypes seen in the samples with Illumina sequence data.

denis_possible_inversions.txt was created as a tab-delimted text file with a line for each putative inversion containing a unique inversion name, chromosome, start, and end. 10 of the random "outlier" windows were also selcted and included to compare against.

#### Map long reads/contigs to reference genome
The steps below can be run on either contigs or reads using an input parameter unless otherwise noted. Many of these scripts take a single sample id as input. This is because the genome assemblies were finishing days apart, and some did not have enough coverage to proudce usable assembled genomes. So it was just easier to run the downstream steps individually for each sample rather than handle filtering and array submissions.

  1. Map read/contig to ref genome using minimap using the inversion_check_minimap.sbatch
  2. Identify the reads/contigs that overlap with putative inversions, and write the locations to a text file that can be subsequently parsed
      * find_inversion_overlaps_in_minimap_results.py
        * Inputs: 
            * sample id
            * sequence type (trimmedReads or contigs)
        
  3. Process the minimap inversion overlap results and store them in a sqlite database (minimap_results.db) 
      * process_inversion_overlaps.py
        * Inputs: 
            * sequence type (trimmedReads or contigs)

  4. Plot the reads/contigs that overlap with putative inversions with inversion_check_contig_map.R. The plot contains two panels:
      * Inputs: 
          * window size
          * sample id
          * sequence type (trimmedReads or contigs)
      * Scatterplot of sech derived allele frequencies where banding of allele frequencies indicates where the putative inversions likely lie.
      * contigs/reads with one per line. Those that map in more than one piece are colored red to draw attention to them.


#### Look for differentiated haplotypes in reads/contigs
Since some of the samples that were sequenced with Oxford Nanopore were not part of the original dataset, we also wanted to check if they had the differentiated haplotypes we had previously seen. Minimap does not return sequences, so LastZ was used to map just the reads/contigs that overlapped with putative inversions. Lastz alignments were then parsed to identify the alleles for sites present in the Illumina sequenced vcf analysis.

Note that this analysis is still ongoing as we await additional Oxford Nanopore sequencing results.

  1. Map the reads/contigs that Minimap found overlapped with putative inversions against the reference genome with lastZ with inversion_check_lastz_reads.sbatch
      * Calls inversion_check_lastz.py
      * Inputs: 
          * samople id
  2. Identify sech anro alleles for polymorphic sites within putative inversions by looking at illumina data with get_inversion_sech_alleles.py
  3. Get the number of derived and sech anro alleles for reads/contigs overlapping putative inversions
      * python3 get_derived_alleles_from_mafs_reads.py 
          * No Inputs 
             
      * python3 get_derived_alleles_from_contig_mafs.py 
          * No Inputs 

  4. Look at anro allele concordance for reads overlapping putative inversions using inversion_check_per_der_reads.R. The R script makes a pdf with plots looking at the histograms of sech anro alleles concordance across all putative inversions and for each putative inversion. Anro allele concordance is just the proportion of alleles that are the same with the representative sech anro allele. The expectation is that a diverged haplotype will have a concordance near zero, and an individuyal heterozygous for the two haplotypes will have an anro allele concordance near 0.5. 
      * No Inputs 

  5. Look at anro allele concordance for reads overlapping putative inversions, coloring reads by their amount of anro allele concordance using inversion_check_per_der_reads_inversion_plots.R. The R script creates a pdf file with one page per putative inversion containing:
      * No Inputs 
      * Scatterplot of sech derived allele frequencies where banding of allele frequencies indicates where the putative inversions likely lie.
      * Scatterplot of the read derived allele frequency. A blue line denotes the average frequency for 2kb windows across the region.
      * Scatterplot of the read sech anro allele frequency. A red line indicates the average frequency for 2kb windows across the region.
      * The mapped reads plotted one per line. The reads are colored by their level of anro allele concordance with blue denoting low concordance, green indicating around 50% concordance, and red denoting high concordance. Extended regions of blue reads should indicate a haplotype diverged from sech anro.



