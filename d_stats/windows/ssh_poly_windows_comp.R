library("RSQLite")  
library(colorRamps)


args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]



#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('ssh_poly_win_comp_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)

#load window data from db
d_stats_table_ssh <- paste('d_stat_win_', win_size, sep='')
d_stats_table_sech_praslin <- paste('d_stat_win_', win_size, '_praslin_vs_sech_base', sep='')
d_stats_table_sech_denis <- paste('d_stat_win_', win_size, '_denis_vs_sech_base', sep='')

poly_table_sech_all <- paste('poly_win_', win_size, sep='')
poly_table_sech_base <- paste('poly_win_', win_size, '_praslin_vs_sech_base', sep='')
poly_table_sech_praslin <- paste('poly_win_', win_size, '_sech_base_vs_praslin', sep='')
poly_table_sech_denis <- paste('poly_win_', win_size, '_sech_base_vs_denis', sep='')

num_sites_cutoff <- as.numeric(win_size) / 50
ab_cutoff <- as.numeric(win_size) / 10000

d_wins_ssh <- dbGetQuery(conn, paste("select *, abba + baba + baaa + abaa ab
                                      from ", d_stats_table_ssh, "
                                      where d_plus is not null
                                      and num_sites > ", num_sites_cutoff, sep=''))

d_wins_sech_praslin <- dbGetQuery(conn, paste("select *, abba + baba + baaa + abaa ab
                                               from ", d_stats_table_sech_praslin, "
                                               where d_plus is not null
                                               and num_sites > ", num_sites_cutoff, sep=''))




poly_wins_sech_all <- dbGetQuery(conn, paste("select * 
                                              from ", poly_table_sech_all, "
                                              where pi_sim is not null
                                              and pi_sech is not null
                                              and pi_ssh is not null
                                              and num_sites > ", num_sites_cutoff, sep=''))

poly_wins_sech_base <- dbGetQuery(conn, paste("select * 
                                               from ", poly_table_sech_base, "
                                               where pi_sim is not null
                                               and pi_sech is not null
                                               and pi_ssh is not null
                                               and num_sites > ", num_sites_cutoff, sep=''))

poly_wins_sech_praslin <- dbGetQuery(conn, paste("select * 
                                                  from ", poly_table_sech_praslin, "
                                                  where pi_sim is not null
                                                  and pi_sech is not null
                                                  and pi_ssh is not null
                                                  and num_sites > ", num_sites_cutoff, sep=''))

poly_wins_sech_denis <- dbGetQuery(conn, paste("select * 
                                                from ", poly_table_sech_denis, "
                                                where pi_sim is not null
                                                and pi_sech is not null
                                                and pi_ssh is not null
                                                and num_sites > ", num_sites_cutoff, sep=''))




chroms <- sort(unique(d_wins_ssh$chrom))
sech_all_threshold <- quantile(poly_wins_sech_all$pi_sech, .99)
pi_range <- range(poly_wins_sech_all$pi_sech, poly_wins_sech_base$pi_sech, poly_wins_sech_praslin$pi_sech, poly_wins_sech_denis$pi_sech)

pdf(pdf_file, height=8, width=10.5)
for (chrom in chroms)
  {
  sech_pi_outlier_filter <- poly_wins_sech_all$chrom == chrom & poly_wins_sech_all$pi_sech > sech_all_threshold
  xlim_range <- c(0, max(poly_wins_sech_all$end[poly_wins_sech_all$chrom == chrom]))  
  

  par(mfrow=c(5,1), mar=c(0, 4.1, 4.1, 2.1))
  win_filter <- d_wins_ssh$chrom == chrom
  ab_range <- range(c(d_wins_ssh$abba,  d_wins_ssh$baba, d_wins_ssh$baaa, d_wins_ssh$abaa))
  plot((d_wins_ssh$start[win_filter] + d_wins_ssh$end[win_filter]) / 2, d_wins_ssh$abba[win_filter], col='darkblue', type='l', ylim=ab_range, xlim=xlim_range, xlab='', xaxt='n', ylab='AB: sech into ssh', main=paste(chrom, "; win size = ", win_size, sep=''))
  points((d_wins_ssh$start[win_filter] + d_wins_ssh$end[win_filter]) / 2, d_wins_ssh$baba[win_filter], col='lightblue', type='l')
  points((d_wins_ssh$start[win_filter] + d_wins_ssh$end[win_filter]) / 2, d_wins_ssh$baaa[win_filter], col='purple', type='l')
  points((d_wins_ssh$start[win_filter] + d_wins_ssh$end[win_filter]) / 2, d_wins_ssh$abaa[win_filter], col='pink', type='l')
  legend('topleft', c('ABBA', 'BABA', 'BAAA', 'ABAA'), fill = c('darkblue', 'lightblue', 'purple', 'pink'), border=NA)


  par(mar=c(0, 4.1, 0, 2.1))
  ab_range <- range(c(d_wins_sech_praslin$abba,  d_wins_sech_praslin$baba, d_wins_sech_praslin$baaa, d_wins_sech_praslin$abaa))
  plot((d_wins_sech_praslin$start[win_filter] + d_wins_sech_praslin$end[win_filter]) / 2, d_wins_sech_praslin$abba[win_filter], col='darkblue', type='l', ylim=ab_range, xlim=xlim_range, xlab='', xaxt='n', ylab='AB: sim into sech_base', main='')
  points((d_wins_sech_praslin$start[win_filter] + d_wins_sech_praslin$end[win_filter]) / 2, d_wins_sech_praslin$baba[win_filter], col='lightblue', type='l')
  points((d_wins_sech_praslin$start[win_filter] + d_wins_sech_praslin$end[win_filter]) / 2, d_wins_sech_praslin$baaa[win_filter], col='purple', type='l')
  points((d_wins_sech_praslin$start[win_filter] + d_wins_sech_praslin$end[win_filter]) / 2, d_wins_sech_praslin$abaa[win_filter], col='pink', type='l')
  legend('topleft', c('ABBA', 'BABA', 'BAAA', 'ABAA'), fill = c('darkblue', 'lightblue', 'purple', 'pink'), border=NA)



  win_filter <- d_wins_ssh$chrom == chrom & d_wins_ssh$ab > ab_cutoff
  plot((d_wins_ssh$start[win_filter] + d_wins_ssh$end[win_filter]) / 2, d_wins_ssh$d_plus[win_filter], type='l', ylim=c(min(d_wins_ssh$d_plus), max(d_wins_ssh$d_plus)), xlim=xlim_range, xlab='', xaxt='n', ylab='D+ sech into ssh', main='')
  abline(h=0, col='grey')
  abline(h=quantile(d_wins_ssh$d_plus, .99), col=adjustcolor('red', .2))
  rect(poly_wins_sech_all$start[sech_pi_outlier_filter], rep(-99, sum(sech_pi_outlier_filter)), poly_wins_sech_all$end[sech_pi_outlier_filter], rep(99, sum(sech_pi_outlier_filter)), col=adjustcolor('darkgreen', .4), border=NA)


  win_filter <- d_wins_sech_praslin$chrom == chrom  & d_wins_sech_praslin$ab > ab_cutoff
  plot((d_wins_sech_praslin$start[win_filter] + d_wins_sech_praslin$end[win_filter]) / 2, d_wins_sech_praslin$d_plus[win_filter], type='l', ylim=c(min(d_wins_sech_praslin$d_plus), max(d_wins_sech_praslin$d_plus)), xlim=xlim_range, xlab='', xaxt='n', ylab='D+ sim into sech_base', main='')
  abline(h=0, col='grey')
  abline(h=quantile(d_wins_sech_praslin$d_plus, .99), col=adjustcolor('red', .2))
  rect(poly_wins_sech_all$start[sech_pi_outlier_filter], rep(-99, sum(sech_pi_outlier_filter)), poly_wins_sech_all$end[sech_pi_outlier_filter], rep(99, sum(sech_pi_outlier_filter)), col=adjustcolor('darkgreen', .4), border=NA)



  win_filter <- poly_wins_sech_all$chrom == chrom
  
  par(mar=c(5.1, 4.1, 0, 2.1))
  plot((poly_wins_sech_all$start[win_filter] + poly_wins_sech_all$end[win_filter]) / 2, poly_wins_sech_all$pi_sech[win_filter], xlim=xlim_range, type='l', col='red', ylim=pi_range, xlab='', ylab='Sech pi', main='')
  points((poly_wins_sech_base$start[win_filter] + poly_wins_sech_base$end[win_filter]) / 2, poly_wins_sech_base$pi_sech[win_filter], type='l', col='orange')
  points((poly_wins_sech_denis$start[win_filter] + poly_wins_sech_denis$end[win_filter]) / 2, poly_wins_sech_denis$pi_sech[win_filter], type='l', col='purple')
  points((poly_wins_sech_praslin$start[win_filter] + poly_wins_sech_praslin$end[win_filter]) / 2, poly_wins_sech_praslin$pi_sech[win_filter], type='l', col='lightblue')
  legend('topleft', c('All', 'Base', 'Denis', 'Praslin'), fill = c('red', 'orange', 'purple', 'lightblue'), border=NA)
  abline(h=sech_all_threshold, col=adjustcolor('darkgreen', .4))

  }
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#test overlap between peaks
boolean_d_plus_outlier <- (d_wins_ssh$d_plus > quantile(d_wins_ssh$d_plus, .99) & d_wins_ssh$ab > ab_cutoff) | (d_wins_sech_praslin$d_plus > quantile(d_wins_sech_praslin$d_plus, .99) & d_wins_sech_praslin$ab > ab_cutoff)
boolean_sech_pi_outlier <- poly_wins_sech_all$pi_sech > sech_all_threshold & d_wins_ssh$ab > ab_cutoff  & d_wins_sech_praslin$ab > ab_cutoff

num_d_plus_outlier <- sum(boolean_d_plus_outlier)
num_sech_pi_outlier <- sum(boolean_sech_pi_outlier)
num_overlap <- sum(boolean_d_plus_outlier & boolean_sech_pi_outlier)
num_win <- sum(d_wins_ssh$ab > ab_cutoff  & d_wins_sech_praslin$ab > ab_cutoff)

paste(num_d_plus_outlier, num_sech_pi_outlier, num_overlap)
#"38 23 10"
dhyper(num_overlap, num_d_plus_outlier, num_win - num_d_plus_outlier, num_sech_pi_outlier)
#[1] 6.977748e-13
