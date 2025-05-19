library("RSQLite")  
library(colorRamps)

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]
table_suffix <- args[2]


#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
gwas_snp_db_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db'
pdf_file <- paste('poly_win_', win_size, '_', table_suffix, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, paste("attach database '", gwas_snp_db_file, "' as g", sep=''))

#define db table
d_stats_table <- paste('d_stat_win_', win_size, '_', table_suffix, sep='')
poly_table <- paste('poly_win_', win_size, sep='')
poly_diff_table <- paste('poly_diff_win_', win_size, sep='')
sfs_table <- paste('sfs_der_freq_', win_size, sep='')

#process pops
pops <- sort(strsplit(table_suffix, '_')[[1]][1:3])
pop_str <- paste("'", paste(pops, collapse = "', '"), "'", sep='')


#load window and gwas data from db
d_wins <- dbGetQuery(conn, paste("select *, abba + baba + baaa + abaa ab
                                  from ", d_stats_table, "
                                  where d_plus is not null
                                  order by chrom, start", sep=''))

poly_wins <- dbGetQuery(conn, paste("select * 
                                     from ", poly_table, " p
                                     where pop in (", pop_str, ")
                                     and exists (select 'x' from ", d_stats_table, " d
                                                 where p.chrom = d.chrom
                                                 and p.start = d.start
                                                 and d_plus is not null)
                                     order by chrom, start, pop", sep=''))


poly_diff_wins <- dbGetQuery(conn, paste("select * 
                                          from ", poly_diff_table, " p
                                          where pop_a in (", pop_str, ")
                                          and pop_b in (", pop_str, ")
                                          and exists (select 'x' from ", d_stats_table, " d
                                                      where p.chrom = d.chrom
                                                      and p.start = d.start
                                                      and d_plus is not null)
                                          order by chrom, start, pop_a, pop_b", sep=''))

#get pop colors
pop_cols <- dbGetQuery(conn, paste("select pop, col
                                       from pop_cols
                                       where pop in (", pop_str, ")
                                       order by pop", sep=''))

col_12 <- colorRampPalette(colors=c(pop_cols$col[1], pop_cols$col[2]))(101)[51]
col_13 <- colorRampPalette(colors=c(pop_cols$col[1], pop_cols$col[3]))(101)[51]
col_23 <- colorRampPalette(colors=c(pop_cols$col[2], pop_cols$col[3]))(101)[51]


#get gwas sites
gwas_pos <- dbGetQuery(conn, paste("select g.chrom, pos, p, start, end, d_plus
                                    from gwas_snp g, ", d_stats_table, " w
                                    where g.chrom = w.chrom
                                    and pos between start and end
                                    and p < 1e-7", sep=''))

#get derived allele site frequency spectrum
der_freq_sfs <- dbGetQuery(conn, paste("select pop, der_freq, num_sites
                                        from ", sfs_table, "
                                        where der_freq > 0
                                        and der_freq < 1
                                        and pop in (", pop_str, ")
                                        order by pop, der_freq", sep=''))

#get unique chromosomes
chroms <- sort(unique(d_wins$chrom))

#get threshold for number of called sites per window
num_sites_cutoff <- as.numeric(win_size) / 5

#get threshold for AB sum per window
ab_cutoff <- as.numeric(win_size) / 10000


pdf(pdf_file, height=8, width=10.5)
#hist of sites per window
hist(d_wins$num_sites, breaks=seq(0, max(d_wins$num_sites) + 100, 100), col='black', xlab='num sites', ylab='windows', main=c('Total sites per window', paste('window size =', win_size)))
abline(v=as.numeric(win_size) / 5, col='red')

par(mfrow=c(3,1))
#hist of pi for all three pops
pi_breaks = seq(min(poly_wins$pi[poly_wins$num_sites > as.numeric(win_size) / 5]), max(poly_wins$pi[poly_wins$num_sites > as.numeric(win_size) / 5]) + .001, .001)
for (i in 1:length(pops))
  {
  hist(poly_wins$pi[poly_wins$num_sites > as.numeric(win_size) / 5 & poly_wins$pop == pops[i]], breaks=pi_breaks, col='black', xlab='pi', ylab='windows', main=c(paste('Pi ', pops[i],  ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size))) 
  }


#hist of derived allele freq for all three pops
der_freq_breaks = seq(min(poly_wins$der_freq[poly_wins$num_sites > as.numeric(win_size) / 5], na.rm = T), max(poly_wins$der_freq[poly_wins$num_sites > as.numeric(win_size) / 5], na.rm = T) + .001, .001)
for (i in 1:length(pops))
  {
  hist(poly_wins$der_freq[poly_wins$num_sites > as.numeric(win_size) / 5 & poly_wins$pop == pops[i]], breaks=der_freq_breaks, col='black', xlab='Derived allele freq', ylab='windows', main=c(paste('Derived freq ', pops[i], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
  }



#define filters for poly_diff_wins dataframe
poly_diff_filter <- poly_diff_wins$num_sites > as.numeric(win_size) / 5
poly_diff_filter_12 <- poly_diff_filter & poly_diff_wins$pop_a == pops[1] & poly_diff_wins$pop_b == pops[2]
poly_diff_filter_13 <- poly_diff_filter & poly_diff_wins$pop_a == pops[1] & poly_diff_wins$pop_b == pops[3]
poly_diff_filter_23 <- poly_diff_filter & poly_diff_wins$pop_a == pops[2] & poly_diff_wins$pop_b == pops[3]

#derived allele freq diff histograms
der_freq_diff_breaks = seq(min(poly_diff_wins$der_freq_diff[poly_diff_filter], na.rm=T), max(poly_diff_wins$der_freq_diff[poly_diff_filter], na.rm=T) + .001, .001)
hist(poly_diff_wins$der_freq_diff[poly_diff_filter_12], breaks=der_freq_diff_breaks, col='black', xlab='Derived allele freq diff', ylab='windows', main=c(paste('Derived freq diff: ', pops[1], ' - ', pops[2], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_diff_wins$der_freq_diff[poly_diff_filter_13], breaks=der_freq_diff_breaks, col='black', xlab='Derived allele freq diff', ylab='windows', main=c(paste('Derived freq diff: ', pops[1], ' - ', pops[3], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_diff_wins$der_freq_diff[poly_diff_filter_23], breaks=der_freq_diff_breaks, col='black', xlab='Derived allele freq diff', ylab='windows', main=c(paste('Derived freq diff: ', pops[2], ' - ', pops[3], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))


#dxy histograms
dxy_breaks = seq(min(poly_diff_wins$dxy[poly_diff_filter]), max(poly_diff_wins$dxy[poly_diff_filter]) + .005, .005)
hist(poly_diff_wins$dxy[poly_diff_filter_12], breaks=dxy_breaks, col='black', xlab='dxy', ylab='windows', main=c(paste('Dxy ', pops[1], ' - ', pops[2], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_diff_wins$dxy[poly_diff_filter_13], breaks=dxy_breaks, col='black', xlab='dxy', ylab='windows', main=c(paste('Dxy ', pops[1], ' - ', pops[3], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_diff_wins$dxy[poly_diff_filter_23], breaks=dxy_breaks, col='black', xlab='dxy', ylab='windows', main=c(paste('Dxy ', pops[2], ' - ', pops[3], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))


#fst histograms
fst_breaks = seq(min(poly_diff_wins$fst[poly_diff_filter]), max(poly_diff_wins$fst[poly_diff_filter]) + .005, .005)
hist(poly_diff_wins$fst[poly_diff_filter_12], breaks=fst_breaks, col='black', xlab='Fst', ylab='windows', main=c(paste('Fst ', pops[1], ' - ', pops[2], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_diff_wins$fst[poly_diff_filter_13], breaks=fst_breaks, col='black', xlab='Fst', ylab='windows', main=c(paste('Fst ', pops[1], ' - ', pops[3], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_diff_wins$fst[poly_diff_filter_23], breaks=fst_breaks, col='black', xlab='Fst', ylab='windows', main=c(paste('Fst ', pops[2], ' - ', pops[3], ' for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))
par(mfrow=c(1,1))



#derived allele site frequency spectrum for all three pops
plot(der_freq_sfs$der_freq[der_freq_sfs$pop == pops[1]], der_freq_sfs$num_sites[der_freq_sfs$pop == pops[1]], type='l', col=pop_cols$col[1], ylim=c(0, max(der_freq_sfs$num_sites[der_freq_sfs$pop %in% pops])), xlab='Derived freq', ylab='Sites', main='Per site derived allele freqs')
points(der_freq_sfs$der_freq[der_freq_sfs$pop == pops[2]], der_freq_sfs$num_sites[der_freq_sfs$pop == pops[2]], type='l', col=pop_cols$col[2])
points(der_freq_sfs$der_freq[der_freq_sfs$pop == pops[3]], der_freq_sfs$num_sites[der_freq_sfs$pop == pops[3]], type='l', col=pop_cols$col[3])
legend('topright', pops, fill = pop_cols$col, border=NA)


#traces of D+ and polymorphism data by windows for each chrom
for (chrom in chroms)
  {
  xlim_range <- c(0, max(poly_wins$end[poly_wins$chrom == chrom]))  
  
  #gwas sites
  par(mfrow=c(7,1), mar=c(0, 4.1, 4.1, 2.1))
  gwas_p_range <- range(-log(gwas_pos$p, 10))
  plot(gwas_pos$pos[gwas_pos$chrom == chrom], -log(gwas_pos$p[gwas_pos$chrom == chrom], 10), pch=20, ylim=gwas_p_range, xlim=xlim_range, xlab='', xaxt='n', ylab='GWAS -log(Pi)', main=paste(chrom, "; win size = ", win_size, sep=''))
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)

  #AB site counts
  par(mar=c(0, 4.1, 0, 2.1))
  win_filter <- d_wins$chrom == chrom & d_wins$num_sites > as.numeric(win_size) / 5
  
  ab_range <- range(c(d_wins$abba[win_filter],  d_wins$baba[win_filter], d_wins$baaa[win_filter], d_wins$abaa[win_filter]))
  plot((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$abba[win_filter], col='darkblue', type='l', ylim=ab_range, xlim=xlim_range, xlab='', xaxt='n', ylab='AB site patterns', main='')
  points((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$baba[win_filter], col='lightblue', type='l')
  points((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$baaa[win_filter], col='purple', type='l')
  points((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$abaa[win_filter], col='pink', type='l')
  legend('topleft', c('ABBA', 'BABA', 'BAAA', 'ABAA'), fill = c('darkblue', 'lightblue', 'purple', 'pink'), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-9999, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(9999, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)


  #D+
  win_filter_d_plus <- d_wins$chrom == chrom & d_wins$num_sites > as.numeric(win_size) / 5 & d_wins$ab > ab_cutoff
  plot((d_wins$start[win_filter_d_plus] + d_wins$end[win_filter_d_plus]) / 2, d_wins$d_plus[win_filter_d_plus], type='l', ylim=c(min(d_wins$d_plus[d_wins$num_sites > as.numeric(win_size) / 5  & d_wins$ab > ab_cutoff]), max(d_wins$d_plus[d_wins$num_sites > as.numeric(win_size) / 5 & d_wins$ab > ab_cutoff])), xlim=xlim_range, xlab='', xaxt='n', ylab='D+', main='')
  abline(h=0, col='grey')
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  abline(h=quantile(d_wins$d_plus[d_wins$num_sites > as.numeric(win_size) / 5], .99), col=adjustcolor('red', .2))

  #define filters for poly_wins dataframe
  win_poly_filter <- poly_wins$chrom == chrom & poly_wins$num_sites > as.numeric(win_size) / 5
  win_poly_filter_1 <- win_poly_filter & poly_wins$pop == pops[1]
  win_poly_filter_2 <- win_poly_filter & poly_wins$pop == pops[2]
  win_poly_filter_3 <- win_poly_filter & poly_wins$pop == pops[3]

  #pi for all three pops
  pi_range <- range(poly_wins$pi[win_poly_filter])
  plot((poly_wins$start[win_poly_filter_1] + poly_wins$end[win_poly_filter_1]) / 2, poly_wins$pi[win_poly_filter_1], xlim=xlim_range, type='l', col=pop_cols$col[1], ylim=pi_range, xlab='', xaxt='n', ylab='Pi', main='')
  points((poly_wins$start[win_poly_filter_2] + poly_wins$end[win_poly_filter_2]) / 2, poly_wins$pi[win_poly_filter_2], type='l', col=pop_cols$col[2])
  points((poly_wins$start[win_poly_filter_3] + poly_wins$end[win_poly_filter_3]) / 2, poly_wins$pi[win_poly_filter_3], type='l', col=pop_cols$col[3])
  legend('topleft', pops, fill = pop_cols$col, border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)


  #define filters for poly_diff_wins dataframe
  win_poly_diff_filter <- poly_diff_wins$num_sites > as.numeric(win_size) / 5 & poly_diff_wins$chrom == chrom
  win_poly_diff_filter_12 <- win_poly_diff_filter & poly_diff_wins$pop_a == pops[1] & poly_diff_wins$pop_b == pops[2]
  win_poly_diff_filter_13 <- win_poly_diff_filter & poly_diff_wins$pop_a == pops[1] & poly_diff_wins$pop_b == pops[3]
  win_poly_diff_filter_23 <- win_poly_diff_filter & poly_diff_wins$pop_a == pops[2] & poly_diff_wins$pop_b == pops[3]
  
  #dxy
  dxy_range <- range(poly_diff_wins$dxy[win_poly_diff_filter])
  plot((poly_diff_wins$start[win_poly_diff_filter_12] + poly_diff_wins$end[win_poly_diff_filter_12]) / 2, poly_diff_wins$dxy[win_poly_diff_filter_12], xlim=xlim_range, type='l', col=col_12, ylim=dxy_range, xlab='', xaxt='n', ylab='Dxy', main='')
  points((poly_diff_wins$start[win_poly_diff_filter_13] + poly_diff_wins$end[win_poly_diff_filter_13]) / 2, poly_diff_wins$dxy[win_poly_diff_filter_13], type='l', col=col_13)
  points((poly_diff_wins$start[win_poly_diff_filter_23] + poly_diff_wins$end[win_poly_diff_filter_23]) / 2, poly_diff_wins$dxy[win_poly_diff_filter_23], type='l', col=col_23)
  legend('topleft', c(paste(pops[1], '-', pops[2], sep=''), paste(pops[1], '-', pops[3], sep=''), paste(pops[2], '-', pops[3], sep='')), fill = c(col_12, col_13, col_23), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  
  #fst
  fst_range <- range(poly_diff_wins$fst[win_poly_diff_filter])
  plot((poly_diff_wins$start[win_poly_diff_filter_12] + poly_diff_wins$end[win_poly_diff_filter_12]) / 2, poly_diff_wins$fst[win_poly_diff_filter_12], xlim=xlim_range, type='l', col=col_12, ylim=fst_range, xlab='Pos', xaxt='n', ylab='Fst', main='')
  points((poly_diff_wins$start[win_poly_diff_filter_13] + poly_diff_wins$end[win_poly_diff_filter_13]) / 2, poly_diff_wins$fst[win_poly_diff_filter_13], type='l', col=col_13)
  points((poly_diff_wins$start[win_poly_diff_filter_23] + poly_diff_wins$end[win_poly_diff_filter_23]) / 2, poly_diff_wins$fst[win_poly_diff_filter_23], type='l', col=col_23)
  legend('topleft', c(paste(pops[1], '-', pops[2], sep=''), paste(pops[1], '-', pops[3], sep=''), paste(pops[2], '-', pops[3], sep='')), fill = c(col_12, col_13, col_23), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  
  #derived allele frequencies for all three populations
  par(mar=c(5.1, 4.1, 0, 2.1))
  der_freq_range <- range(poly_wins$der_freq[win_poly_filter], na.rm=T)
  plot((poly_wins$start[win_poly_filter_1] + poly_wins$end[win_poly_filter_1]) / 2, poly_wins$der_freq[win_poly_filter_1], xlim=xlim_range, type='l', col='blue', ylim=der_freq_range, xlab='', ylab='Der freq', main='')
  points((poly_wins$start[win_poly_filter_2] + poly_wins$end[win_poly_filter_2]) / 2, poly_wins$der_freq[win_poly_filter_2], type='l', col='red')
  points((poly_wins$start[win_poly_filter_3] + poly_wins$end[win_poly_filter_3]) / 2, poly_wins$der_freq[win_poly_filter_3], type='l', col='gold')
  legend('topleft', pops, fill = pop_cols$col, border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  }
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))


#pi sech specific plots to help better understand how pi sech relates to introgression stats between sim and sech
if ('sech' %in% pops)
  {
  win_filter_d_plus <- d_wins$num_sites > as.numeric(win_size) / 5
  poly_filter <- poly_wins$num_sites > as.numeric(win_size) / 5
  
  #scatter plot of pi sech vs D+ (each point is a window)
  plot(poly_wins$pi[poly_filter & poly_wins$pop == 'sech'], d_wins$d_plus[win_filter_d_plus], pch=20, col=adjustcolor('gray', .6), xlab='pi sech', ylab='D+', main='pi sech vs D+')

  #scatter plot of pi sech vs BAAA - ABAA with points colored by D+ value
  d_plus_col_pal <- adjustcolor(matlab.like(101), .4)
  d_plus_cols <- d_plus_col_pal[1 + round(100 * (d_wins$d_plus[win_filter_d_plus] - min(d_wins$d_plus[win_filter_d_plus])) / (max(d_wins$d_plus[win_filter_d_plus]) - min(d_wins$d_plus[win_filter_d_plus])))]
  plot(poly_wins$pi[poly_filter & poly_wins$pop == 'sech'], d_wins$baaa[win_filter_d_plus] - d_wins$abaa[win_filter_d_plus], pch=20, col=d_plus_cols, xlab='pi sech', ylab='BAAA-ABAA', main='pi sech vs BAAA - ABAA')
  
  #plot custom legend
  x_len <- max(poly_wins$pi[poly_filter & poly_wins$pop == 'sech']) - min(poly_wins$pi[poly_filter & poly_wins$pop == 'sech']) 
  y_len <- max((d_wins$baaa - d_wins$abaa)[win_filter_d_plus]) - min((d_wins$baaa - d_wins$abaa)[win_filter_d_plus]) 
  
  x_min <- min(poly_wins$pi[poly_filter & poly_wins$pop == 'sech']) + .9 * x_len
  x_max <- x_min + 0.05 * x_len
  y_min <- min((d_wins$baaa - d_wins$abaa)[win_filter_d_plus]) + 0.1 * y_len
  y_max <- y_min + 0.25 * y_len
  y_step <- (y_max - y_min) / length(d_plus_col_pal)
  
  for (i in 1:length(d_plus_cols))
    {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=d_plus_col_pal[i], border=d_plus_col_pal[i])}
    
  legend_d_plus <- round(seq(min(round(d_wins$d_plus[win_filter_d_plus], 2)), max(round(d_wins$d_plus[win_filter_d_plus], 2)), length.out=5),2)
  text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_d_plus)), legend_d_plus, cex=.6)
  text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'D+')
  
  
  #scatter plot of D+ vs BAAA-ABAA with points colored by pi sech and point size scales with pi sech to help outliers standout
  sech_pi_col_pal <- adjustcolor(matlab.like(101), .4)
  sech_pi_cols <- sech_pi_col_pal[1 + round(100 * (poly_wins$pi[poly_filter & poly_wins$pop == 'sech'] - min(poly_wins$pi[poly_filter & poly_wins$pop == 'sech'])) / (max(poly_wins$pi[poly_filter & poly_wins$pop == 'sech']) - min(poly_wins$pi[poly_filter & poly_wins$pop == 'sech'])))]
  cexes <- 3 * poly_wins$pi[poly_filter & poly_wins$pop == 'sech'] / max(poly_wins$pi[poly_filter & poly_wins$pop == 'sech'])
  plot(d_wins$baaa[win_filter_d_plus] - d_wins$abaa[win_filter_d_plus], d_wins$d_plus[win_filter_d_plus], pch=20, cex=cexes, col=sech_pi_cols, xlab='BAAA-ABAA', ylab='D+', main='BAAA - ABAA vs D+')
  
  #plot custom legend
  x_len <- max((d_wins$baaa - d_wins$abaa)[win_filter_d_plus]) - min((d_wins$baaa - d_wins$abaa)[win_filter_d_plus]) 
  y_len <- max(d_wins$d_plus[win_filter_d_plus]) - min(d_wins$d_plus[win_filter_d_plus]) 
  
  x_min <- min((d_wins$baaa - d_wins$abaa)[win_filter_d_plus]) + 0.9 * x_len
  x_max <- x_min + 0.05 * x_len
  y_min <- min(d_wins$d_plus[win_filter_d_plus]) + .05 * y_len
  y_max <- y_min + 0.25 * y_len
  y_step <- (y_max - y_min) / length(sech_pi_col_pal)
  
  for (i in 1:length(sech_pi_cols))
    {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=sech_pi_col_pal[i], border=sech_pi_col_pal[i])}
    
  legend_pi_sech <- round(seq(min(round(poly_wins$pi[win_filter_d_plus & poly_wins$pop == 'sech'], 2)), max(round(poly_wins$pi[win_filter_d_plus & poly_wins$pop == 'sech'], 2)), length.out=5),2)
  text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_pi_sech)), legend_pi_sech, cex=.6)
  text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'pi sech')
  
  
  #scatter plot of D+ vs ABBA-BABA with points colored by pi sech and point size scales with pi sech to help outliers standout
  cexes <- 3 * poly_wins$pi[win_filter_d_plus & poly_wins$pop == 'sech'] / max(poly_wins$pi[win_filter_d_plus & poly_wins$pop == 'sech'])
  plot(d_wins$abba[win_filter_d_plus] - d_wins$baba[win_filter_d_plus], d_wins$d_plus[win_filter_d_plus], pch=20, cex=cexes, col=sech_pi_cols, xlab='ABBA-BABA', ylab='D+', main='ABBA - BABA vs D+')

  #plot custom legend
  x_len <- max((d_wins$abba - d_wins$baba)[win_filter_d_plus]) - min((d_wins$abba - d_wins$baba)[win_filter_d_plus]) 
  y_len <- max(d_wins$d_plus[win_filter_d_plus]) - min(d_wins$d_plus[win_filter_d_plus]) 
  
  x_min <- min((d_wins$abba - d_wins$baba)[win_filter_d_plus]) + 0.9 * x_len
  x_max <- x_min + 0.05 * x_len
  y_min <- min(d_wins$d_plus[win_filter_d_plus]) + .05 * y_len
  y_max <- y_min + 0.25 * y_len
  y_step <- (y_max - y_min) / length(sech_pi_col_pal)
  
  for (i in 1:length(sech_pi_cols))
    {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=sech_pi_col_pal[i], border=sech_pi_col_pal[i])}
    
  legend_pi_sech <- round(seq(min(round(poly_wins$pi[win_filter_d_plus & poly_wins$pop == 'sech'], 2)), max(round(poly_wins$pi[win_filter_d_plus & poly_wins$pop == 'sech'], 2)), length.out=5),2)
  text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_pi_sech)), legend_pi_sech, cex=.6)
  text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'pi sech')
  }
dev.off()

