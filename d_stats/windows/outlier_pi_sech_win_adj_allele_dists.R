library("RSQLite")  
library(fields)
library(colorRamps)

args <- commandArgs(trailingOnly = TRUE)
win_size <- as.character(args[1])

#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('ssh_outlier_pi_sech_adj_allele_dists_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)

allele_dist_d_plus_table <- paste('outlier_d_plus_win_sech_adj_allele_dist_', win_size, sep='')
allele_dist_pi_sech_table <- paste('outlier_pi_sech_win_sech_adj_allele_dist_', win_size, sep='')
allele_dist_random_table <- paste('outlier_random_win_sech_adj_allele_dist_', win_size, sep='')



allele_dist_d_plus <- dbGetQuery(conn, paste("select *
                                              from ", allele_dist_d_plus_table, sep=''))

allele_dist_pi_sech <- dbGetQuery(conn, paste("select *
                                               from ", allele_dist_pi_sech_table, sep=''))

allele_dist_random <- dbGetQuery(conn, paste("select *
                                              from ", allele_dist_random_table, sep=''))


samples <- dbGetQuery(conn, paste("select location, sample_id
                                   from sample_pop
                                   where pop = 'sechellia'
                                   order by location, sample_id", sep=''))


sech_loc <- data.frame(cbind(c('sech - Anro', 'sech - Denis', 'sech - La Digue', 'sech - Marianne', 'sech - Praslin', 'sech - Unk'),
                 c('Anro, Seychelles', 'Denis, Seychelles', 'La Digue, Seychelles', 'Marianne, Seychelles', 'Praslin, Seychelles', 'Unknown'),
                 c('orange', 'purple', 'green', 'darkblue', 'lightblue', 'gray')))
names(sech_loc) <- c('display_loc', 'location', 'loc_col')


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

