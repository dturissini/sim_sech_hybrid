library(colorRamps)

setwd('/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows')

chroms <- c('2L', '2R', '3L', '3R', 'X')
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'

win_file <- paste(base_dir, 'sim_win.txt', sep='')
pattern_file <- paste(base_dir, 'ssh_win_50k_d_stats.txt', sep='')
gwas_results_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_gwas.oa_resistance.glm.linear'
 
wins <- read.table(win_file, header=T)
patterns <- read.table(pattern_file, header=F, col.names=c('ABBA', 'BABA', 'BBAA', 'BAAA', 'ABAA', 'AABA'))

gwas <- read.table(gwas_results_file, comment.char='', header=T)
names(gwas)[1] <- 'CHROM'
gwas_filtered <- subset(gwas, gwas$P < 1e-7)


d_plus <- (patterns$ABBA - patterns$BABA + patterns$BAAA - patterns$ABAA) / (patterns$ABBA + patterns$BABA + patterns$BAAA + patterns$ABAA)
per_abba <- 100 * patterns$ABBA / (patterns$ABBA + patterns$BABA + patterns$BAAA + patterns$ABAA)
per_baaa <- 100 * patterns$BAAA / (patterns$ABBA + patterns$BABA + patterns$BAAA + patterns$ABAA)
abba_min_baba <- patterns$ABBA - patterns$BABA
baaa_min_abaa <- patterns$BAAA - patterns$ABAA

pdf('ssh_dplus_windows.pdf', height=8, width=10.5)
hist(wins$tot_sites, breaks=seq(0, max(wins$tot_sites) + 100, 100), col='black', xlab='num sites', ylab='windows', main='Total sites per window')
abline(v=1000, col='red')

hist(d_plus[wins$tot_sites > 1000], breaks=100, col='black', xlab='D+', ylab='windows', main='D+ for windows with > 1000 sites')

for (chrom in chroms)
  {
  win_filter <- wins$chr == chrom & wins$tot_sites > 1000
  plot((wins$start[win_filter] + wins$end[win_filter]) / 2, d_plus[win_filter], type='l', ylim=c(min(d_plus[wins$tot_sites > 1000]), max(d_plus[wins$tot_sites > 1000])), xlab='pos', ylab='D+', main=c('D+ 50kb windows', chrom))
  abline(h=0, col='grey')

  chrom_i <- which(c('2L', '2R', '3L', '3R', '4', 'X') == chrom)
  points(gwas_filtered$POS[gwas_filtered$CHROM == chrom_i], rep(max(d_plus[wins$tot_sites > 1000]), sum(gwas_filtered$CHROM == chrom_i)), pch=20, col='red') 
  text(min(wins$start[win_filter]), .95 * max(d_plus[wins$tot_sites > 1000]), 'GWAS', col='red')
  }


plot(d_plus[wins$tot_sites > 1000], abba_min_baba[wins$tot_sites > 1000], pch=20, xlab='D+', ylab='ABBA - BABA', main='D+ vs ABBA - BABA for 50kb windows')
plot(d_plus[wins$tot_sites > 1000], baaa_min_abaa[wins$tot_sites > 1000], pch=20, xlab='D+', ylab='BAAA - ABAA', main='D+ vs BAAA - ABAA for 50kb windows')

d_plus_col_pal <- adjustcolor(matlab.like(101), .4)
d_plus_cols <- d_plus_col_pal[1 + round(100 * (d_plus[wins$tot_sites > 1000] - min(d_plus[wins$tot_sites > 1000])) / (max(d_plus[wins$tot_sites > 1000]) - min(d_plus[wins$tot_sites > 1000])))]
plot(abba_min_baba[wins$tot_sites > 1000], baaa_min_abaa[wins$tot_sites > 1000], pch=20, col=d_plus_cols, xlab='ABBA - BABA', ylab='BAAA - ABAA', main='ABBA - BABA vs BAAA - ABAA for 50kb windows')

abline(v=0, col='gray')
abline(h=0, col='gray')


#plot custom legend
x_len <- max(abba_min_baba[wins$tot_sites > 1000]) - min(abba_min_baba[wins$tot_sites > 1000]) 
y_len <- max(baaa_min_abaa[wins$tot_sites > 1000]) - min(baaa_min_abaa[wins$tot_sites > 1000]) 

x_min <- min(abba_min_baba[wins$tot_sites > 1000])
x_max <- x_min + 0.05 * x_len
y_min <- min(baaa_min_abaa[wins$tot_sites > 1000]) + 0.7 * y_len
y_max <- y_min + 0.25 * y_len

y_step <- (y_max - y_min) / length(d_plus_col_pal)

for (i in 1:length(d_plus_cols))
  {
  rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=d_plus_col_pal[i], border=d_plus_col_pal[i])  
  }
  
legend_d_plus <- round(seq(min(round(d_plus[wins$tot_sites > 1000], 2)), max(round(d_plus[wins$tot_sites > 1000], 2)), length.out=5),2)
text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_d_plus)), legend_d_plus, cex=.6)
text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'D+')



d_plus_binary_cols <- rep(adjustcolor('black', .4), sum(wins$tot_sites > 1000))
d_plus_binary_cols[d_plus[wins$tot_sites > 1000] < 0] <- adjustcolor('red', .4)
plot(abba_min_baba[wins$tot_sites > 1000], baaa_min_abaa[wins$tot_sites > 1000], pch=20, col=d_plus_binary_cols, xlab='ABBA - BABA', ylab='BAAA - ABAA', main='ABBA - BABA vs BAAA - ABAA for 50kb windows')
abline(v=0, col='gray')
abline(h=0, col='gray')
legend('topleft', c('D+ >= 0', 'D+ < 0'), fill=c('black', 'red'), border=NA)


plot(d_plus[wins$tot_sites > 1000], per_abba[wins$tot_sites > 1000], pch=20, xlab='D+', ylab='% ABBA sites', main='D+ vs % ABBA sites for 50kb windows')

for (chrom in chroms)
  {
  win_filter <- wins$chr == chrom & wins$tot_sites > 1000
  plot((wins$start[win_filter] + wins$end[win_filter]) / 2, per_abba[win_filter], type='l', ylim=c(0, 100), xlab='pos', ylab='% ABBA sites', main=c('% ABBA sitres 50kb windows', chrom))
  abline(h=0, col='grey')

  chrom_i <- which(c('2L', '2R', '3L', '3R', '4', 'X') == chrom)
  points(gwas_filtered$POS[gwas_filtered$CHROM == chrom_i], rep(max(d_plus[wins$tot_sites > 1000]), sum(gwas_filtered$CHROM == chrom_i)), pch=20, col='red') 
  text(min(wins$start[win_filter]), .95 * max(d_plus[wins$tot_sites > 1000]), 'GWAS', col='red')
  }

plot(d_plus[wins$tot_sites > 1000], per_baaa[wins$tot_sites > 1000], pch=20, xlab='D+', ylab='% BAAA sites', main='D+ vs % BAAA sites for 50kb windows')

for (chrom in chroms)
  {
  win_filter <- wins$chr == chrom & wins$tot_sites > 1000
  plot((wins$start[win_filter] + wins$end[win_filter]) / 2, per_baaa[win_filter], type='l', ylim=c(0, 100), xlab='pos', ylab='% BAAA sites', main=c('% BAAA sitres 50kb windows', chrom))
  abline(h=0, col='grey')

  chrom_i <- which(c('2L', '2R', '3L', '3R', '4', 'X') == chrom)
  points(gwas_filtered$POS[gwas_filtered$CHROM == chrom_i], rep(max(d_plus[wins$tot_sites > 1000]), sum(gwas_filtered$CHROM == chrom_i)), pch=20, col='red') 
  text(min(wins$start[win_filter]), .95 * max(d_plus[wins$tot_sites > 1000]), 'GWAS', col='red')
  }

dev.off()


which(wins$start < gwas_filtered$POS[gwas_filtered$CHROM == 2 & gwas_filtered$P < 1e-12])