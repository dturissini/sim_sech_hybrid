library("RSQLite")  
library(colorRamps)

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- as.character(args[1])
pop_str <- as.character(args[2])


#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('outlier_pi_sech_adj_allele_dists_', win_size, '_', pop_str, '.pdf', sep='')

setwd(base_dir)


#make db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)


#define db tables
allele_dist_d_plus_table <- paste('outlier_d_plus_win_sech_adj_allele_dist_', win_size, '_', pop_str, sep='')
allele_dist_pi_sech_table <- paste('outlier_pi_sech_win_sech_adj_allele_dist_', win_size, '_', pop_str, sep='')
allele_dist_random_table <- paste('outlier_random_win_sech_adj_allele_dist_', win_size, '_', pop_str, sep='')

#define str of pops to use in db queries
pop_sql_str <- "'sechanro', 'sechden', 'sechlad', 'sechmari', 'sechpras', 'sechunk'"


#get allele distances
allele_dist_d_plus <- dbGetQuery(conn, paste("select *
                                              from ", allele_dist_d_plus_table, sep=''))

allele_dist_pi_sech <- dbGetQuery(conn, paste("select *
                                               from ", allele_dist_pi_sech_table, sep=''))

allele_dist_random <- dbGetQuery(conn, paste("select *
                                              from ", allele_dist_random_table, sep=''))


#get samples and pops
samples <- dbGetQuery(conn, paste("select pop, sample_id
                                   from sample_pop_link
                                   where pop in (", pop_sql_str, ")
                                   order by pop, sample_id", sep=''))



#make plots of the histograms of adjusted allele distances for each sample
pdf(pdf_file, height=8, width=10.5)
for (i in 1:nrow(samples))
  {
  hist_dp <- hist(allele_dist_d_plus$adj_allele_dist[allele_dist_d_plus$sample_id == samples$sample_id[i]], breaks = seq(0, 110, 10), right=F, plot=F)
  hist_ps <- hist(allele_dist_pi_sech$adj_allele_dist[allele_dist_pi_sech$sample_id == samples$sample_id[i]], breaks = seq(0, 110, 10), right=F, plot=F)
  hist_r <- hist(allele_dist_random$adj_allele_dist[allele_dist_random$sample_id == samples$sample_id[i]], breaks = seq(0, 110, 10), right=F, plot=F)

  plot(hist_dp$mids, hist_dp$counts, type='l', col='blue', ylim=c(0, max(c(hist_dp$counts, hist_ps$counts, hist_r$counts))), xlab='adjusted allele distance', ylab='windows', main=c(samples$sample_id[i], samples$location[i]))
  points(hist_ps$mids, hist_ps$counts, type='l', col='red')
  points(hist_r$mids, hist_r$counts, type='l', col='black')
  
  legend('topright', c('D+ outlier', 'pi sech outlier', 'random'), fill=c('blue', 'red', 'black'), border=NA)
  }
dev.off()

