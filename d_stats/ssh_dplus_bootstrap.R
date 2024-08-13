setwd('/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats')

chroms <- c('2L', '2R', '3L', '3R', 'X')
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/'

pdf('ssh_dplus_bootstrap.pdf', height=8, width=16)
par(mfrow=c(1,2))
for (chrom in chroms)
  {
  d_file <- paste(base_dir, 'ripta_', chrom, '/simulans_sim_sech_hybrid_sechellia_melanogaster_observed_values.txt', sep='')
  bootstrap_file <- paste(base_dir, 'ripta_', chrom, '/simulans_sim_sech_hybrid_sechellia_melanogaster_bootstraps.txt', sep='')
  
  d <- read.table(d_file, header=T)
  bootstraps <- read.table(bootstrap_file, header=T)
  
  hist(bootstraps$D, col='black', xlab='D', ylab='bootstraps', main=c('D', chrom)) 
  abline(v=d$D, col='red') 
  
  hist(bootstraps$D., col='black', xlab='D+', ylab='bootstraps', main=c('D+', paste('ABBA - BABA = ', round(d$ABBA - d$BABA,1), ', BAAA - ABAA = ', round(d$BAAA - d$ABAA,1), sep='')))  
  abline(v=d$D., col='red') 
  }
par(mfrow=c(1,1))
dev.off()


