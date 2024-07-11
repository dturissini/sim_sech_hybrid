library("RSQLite")  
library(colorRamps)


args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]
table_suffix_a <- args[2]
table_suffix_b <- args[3]



#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('ssh_poly_win_comp_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)


pops_a <- sort(strsplit(table_suffix_a, '_')[[1]][1:3])

pops_b <- sort(strsplit(table_suffix_b, '_')[[1]][1:3])

pops <- unique(c(pops_a, pops_b))
pop_str <- paste("'", paste(pops, collapse = "', '"), "'", sep='')

sech_pops <- pops[substr(pops, 1, 4) == 'sech']



#load window data from db
d_stats_table_a <- paste('d_stat_win_', win_size, '_', table_suffix_a, sep='')
d_stats_table_b <- paste('d_stat_win_', win_size, '_', table_suffix_b, sep='')


poly_table <- paste('poly_win_', win_size, sep='')


num_sites_cutoff <- as.numeric(win_size) / 50
ab_cutoff <- as.numeric(win_size) / 10000

d_wins_a <- dbGetQuery(conn, paste("select *, abba + baba + baaa + abaa ab
                                      from ", d_stats_table_a, "
                                      where d_plus is not null
                                      and num_sites > ", num_sites_cutoff, sep=''))

d_wins_b <- dbGetQuery(conn, paste("select *, abba + baba + baaa + abaa ab
                                      from ", d_stats_table_b, "
                                      where d_plus is not null
                                      and num_sites > ", num_sites_cutoff, sep=''))




poly_wins <- dbGetQuery(conn, paste("select * 
                                     from ", poly_table, "
                                     where pop in (", pop_str, ")
                                     and num_sites > ", num_sites_cutoff, sep=''))


pop_cols <- dbGetQuery(conn, paste("select pop, col
                                       from pop_cols
                                       where pop in (", pop_str, ")
                                       order by pop", sep=''))


chroms <- sort(unique(d_wins_a$chrom))
sech_all_threshold <- quantile(poly_wins$pi[poly_wins$pop == 'sech'], .99)
pi_range <- range(poly_wins$pi[substr(poly_wins$pop, 1, 4) == 'sech'])

pdf(pdf_file, height=8, width=10.5)
for (chrom in chroms)
  {
  sech_pi_outlier_filter <- poly_wins$chrom == chrom & poly_wins$pi[poly_wins$pop == 'sech'] > sech_all_threshold
  xlim_range <- c(0, max(poly_wins$end[poly_wins$chrom == chrom]))  
  

  par(mfrow=c(5,1), mar=c(0, 4.1, 4.1, 2.1))
  win_filter <- d_wins_a$chrom == chrom
  ab_range <- range(c(d_wins_a$abba,  d_wins_a$baba, d_wins_a$baaa, d_wins_a$abaa))
  plot((d_wins_a$start[win_filter] + d_wins_a$end[win_filter]) / 2, d_wins_a$abba[win_filter], col='darkblue', type='l', ylim=ab_range, xlim=xlim_range, xlab='', xaxt='n', ylab=table_suffix_a, main=paste(chrom, "; win size = ", win_size, sep=''))
  points((d_wins_a$start[win_filter] + d_wins_a$end[win_filter]) / 2, d_wins_a$baba[win_filter], col='lightblue', type='l')
  points((d_wins_a$start[win_filter] + d_wins_a$end[win_filter]) / 2, d_wins_a$baaa[win_filter], col='purple', type='l')
  points((d_wins_a$start[win_filter] + d_wins_a$end[win_filter]) / 2, d_wins_a$abaa[win_filter], col='pink', type='l')
  legend('topleft', c('ABBA', 'BABA', 'BAAA', 'ABAA'), fill = c('darkblue', 'lightblue', 'purple', 'pink'), border=NA)


  par(mar=c(0, 4.1, 0, 2.1))
  ab_range <- range(c(d_wins_b$abba,  d_wins_b$baba, d_wins_b$baaa, d_wins_b$abaa))
  plot((d_wins_b$start[win_filter] + d_wins_b$end[win_filter]) / 2, d_wins_b$abba[win_filter], col='darkblue', type='l', ylim=ab_range, xlim=xlim_range, xlab='', xaxt='n', ylab=table_suffix_b, main='')
  points((d_wins_b$start[win_filter] + d_wins_b$end[win_filter]) / 2, d_wins_b$baba[win_filter], col='lightblue', type='l')
  points((d_wins_b$start[win_filter] + d_wins_b$end[win_filter]) / 2, d_wins_b$baaa[win_filter], col='purple', type='l')
  points((d_wins_b$start[win_filter] + d_wins_b$end[win_filter]) / 2, d_wins_b$abaa[win_filter], col='pink', type='l')
  legend('topleft', c('ABBA', 'BABA', 'BAAA', 'ABAA'), fill = c('darkblue', 'lightblue', 'purple', 'pink'), border=NA)



  win_filter <- d_wins_a$chrom == chrom & d_wins_a$ab > ab_cutoff
  plot((d_wins_a$start[win_filter] + d_wins_a$end[win_filter]) / 2, d_wins_a$d_plus[win_filter], type='l', ylim=c(min(d_wins_a$d_plus), max(d_wins_a$d_plus)), xlim=xlim_range, xlab='', xaxt='n', ylab=table_suffix_a, main='')
  abline(h=0, col='grey')
  abline(h=quantile(d_wins_a$d_plus, .99), col=adjustcolor('red', .2))
  rect(poly_wins$start[sech_pi_outlier_filter], rep(-99, sum(sech_pi_outlier_filter)), poly_wins$end[sech_pi_outlier_filter], rep(99, sum(sech_pi_outlier_filter)), col=adjustcolor('darkgreen', .2), border=NA)


  win_filter <- d_wins_b$chrom == chrom  & d_wins_b$ab > ab_cutoff
  plot((d_wins_b$start[win_filter] + d_wins_b$end[win_filter]) / 2, d_wins_b$d_plus[win_filter], type='l', ylim=c(min(d_wins_b$d_plus), max(d_wins_b$d_plus)), xlim=xlim_range, xlab='', xaxt='n', ylab=table_suffix_b, main='')
  abline(h=0, col='grey')
  abline(h=quantile(d_wins_b$d_plus, .99), col=adjustcolor('red', .2))
  rect(poly_wins$start[sech_pi_outlier_filter], rep(-99, sum(sech_pi_outlier_filter)), poly_wins$end[sech_pi_outlier_filter], rep(99, sum(sech_pi_outlier_filter)), col=adjustcolor('darkgreen', .2), border=NA)



  
  par(mar=c(5.1, 4.1, 0, 2.1))
  plot(1, type='n', xlim=xlim_range, ylim=pi_range, xlab='', ylab='Sech pi', main='')
  for (pop in sech_pops)
    {
    points((poly_wins$start[poly_wins$chrom == chrom & poly_wins$pop == pop] + poly_wins$end[poly_wins$chrom == chrom & poly_wins$pop == pop]) / 2, poly_wins$pi[poly_wins$chrom == chrom & poly_wins$pop == pop], type='l', col=pop_cols$col[pop_cols$pop == pop])
    }
  
  legend('topleft', sech_pops, fill = pop_cols$col[pop_cols$pop %in% sech_pops], border=NA)
  abline(h=sech_all_threshold, col=adjustcolor('darkgreen', .2))

  }
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#test overlap between peaks
quantiles <- c(.99, .98, .95)
for (quantile in quantiles)
  {
  boolean_d_plus_outlier <- (d_wins_a$d_plus > quantile(d_wins_a$d_plus, quantile) & d_wins_a$ab > ab_cutoff) | (d_wins_b$d_plus > quantile(d_wins_b$d_plus, quantile) & d_wins_b$ab > ab_cutoff)
  boolean_sech_pi_outlier <- poly_wins$pop == 'sech' & poly_wins$pi > quantile(poly_wins$pi[poly_wins$pop == 'sech'], quantile) & d_wins_a$ab > ab_cutoff  & d_wins_b$ab > ab_cutoff
  
  num_d_plus_outlier <- sum(boolean_d_plus_outlier)
  num_sech_pi_outlier <- sum(boolean_sech_pi_outlier)
  num_overlap <- sum(boolean_d_plus_outlier & boolean_sech_pi_outlier)
  num_win <- sum(d_wins_a$ab > ab_cutoff  & d_wins_b$ab > ab_cutoff)
  
  cat("quantile", quantile, "\n")
  cat(num_d_plus_outlier, num_sech_pi_outlier, num_overlap, "\n")
  cat(dhyper(num_overlap, num_d_plus_outlier, num_win - num_d_plus_outlier, num_sech_pi_outlier), "\n\n\n")
  }
  
##quantile 0.99 
##38 23 11 
##1.084038e-14 
##
##
##quantile 0.98 
##80 45 19 
##6.524933e-17 
##
##
##quantile 0.95 
##201 111 40 
##7.746315e-16 
