library("RSQLite")  
library(colorRamps)


#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]
pop_str <- args[2]


#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
gwas_snp_db_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db'
pdf_file <- paste('dplus_win_', win_size, '_', pop_str, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, paste("attach database '", gwas_snp_db_file, "' as g", sep=''))


#load window and gwas data from db
d_stats_table <- paste('d_stat_win_', win_size, '_', pop_str, sep='')
wins <- dbGetQuery(conn, paste("select *, abba + baba + baaa + abaa ab, abba + baba abba_baba, baaa + abaa baaa_abaa
                                from ", d_stats_table, "
                                where d_plus is not null", sep=''))

gwas_pos <- dbGetQuery(conn, paste("select g.chrom, pos, p, start, end, d_plus
                                    from gwas_snp g, ", d_stats_table, " w
                                    where g.chrom = w.chrom
                                    and pos between start and end
                                    and p < 1e-7", sep=''))

#get unique chromosomes
chroms <- sort(unique(wins$chrom))

#set threshold for number of called sites per window
as.numeric(win_size) / 5 <- as.numeric(win_size) / 50


#set threshold for ab sum per window
ab_cutoff <- as.numeric(win_size) / 10000


pdf(pdf_file, height=8, width=10.5)
#hist of sites per window
hist(wins$num_sites, breaks=seq(0, max(wins$num_sites) + 100, 100), col='black', xlab='num sites', ylab='windows', main=c('Total sites per window', paste('window size =', win_size)))
abline(v=as.numeric(win_size) / 5, col='red')


#histograms of AB totals with and without cutoffs
par(mfrow=c(3,1))
hist(wins$ab, breaks=seq(0, max(wins$ab) + 10, 10), col='black', xlab='AB total', ylab='windows', main=c('AB total', paste('window size =', win_size)))
abline(v=ab_cutoff, col='red')

hist(wins$abba_baba, breaks=seq(0, max(wins$ab) + 10, 10), col='black', xlab='ABBA_BABA total', ylab='windows', main=c('ABBA_BABA total', paste('window size =', win_size)))
hist(wins$baaa_abaa, breaks=seq(0, max(wins$ab) + 10, 10), col='black', xlab='BAAA_ABAA total', ylab='windows', main=c('BAAA_ABAA total', paste('window size =', win_size)))


hist(wins$ab[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], breaks=seq(0, max(wins$ab) + 10, 10), col='black', xlab='AB total', ylab='windows', main=c(paste('AB total (num_sites > ', as.numeric(win_size) / 5, '; AB total > ', ab_cutoff, ')', sep=''), paste('window size =', win_size)))
abline(v=ab_cutoff, col='red')

hist(wins$abba_baba[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], breaks=seq(0, max(wins$ab) + 10, 10), col='black', xlab='ABBA_BABA total', ylab='windows', main=c(paste('ABBA+BABA total (num_sites > ', as.numeric(win_size) / 5, '; AB total > ', ab_cutoff, ')', sep=''), paste('window size =', win_size)))
hist(wins$baaa_abaa[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], breaks=seq(0, max(wins$ab) + 10, 10), col='black', xlab='BAAA_ABAA total', ylab='windows', main=c(paste('BAAA+ABAA total (num_sites > ', as.numeric(win_size) / 5, '; AB total > ', ab_cutoff, ')', sep=''), paste('window size =', win_size)))
par(mfrow=c(1,1))


#scatterplot of number of called sites and AB sum by window
plot(wins$num_sites, wins$ab, pch=20, col=adjustcolor('gray', .6), xlab='num sites', ylab='AB total', main='num sites vs AB total')
abline(v=as.numeric(win_size) / 5, col='red')
abline(h=ab_cutoff, col='red')



#hist of D+
d_plus_hist <- hist(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5  & wins$ab > ab_cutoff], breaks=100, col='black', xlab='D+', ylab='windows', main=c(paste('D+ for windows with num_sites > ', as.numeric(win_size) / 5, ' and AB total > ', ab_cutoff, sep=''), paste('window size =', win_size)))
abline(v=quantile(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], .99), col='red')
text(quantile(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], .99) + .1, max(d_plus_hist$counts), paste('.99 quantile =', round(quantile(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], .99), 3)), col='red')


#scatterplots of D+ vs AB sum to help identify a relevant ab threshold
par(mfrow=c(3,1))
plot(wins$ab[wins$num_sites > as.numeric(win_size) / 5], wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], pch=20, col=adjustcolor('gray', .6), xlim=c(0, max(wins$ab[wins$num_sites > as.numeric(win_size) / 5])), xlab='AB total', ylab='D+', main=c('AB total vs D+', paste('num sites >', as.numeric(win_size) / 5)))
abline(v=ab_cutoff, col='red')
abline(h=quantile(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], .99), col='red')

plot(wins$abba_baba[wins$num_sites > as.numeric(win_size) / 5], wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], pch=20, col=adjustcolor('gray', .6), xlim=c(0, max(wins$ab[wins$num_sites > as.numeric(win_size) / 5])), xlab='ABBA_BABA total', ylab='D+', main=c('ABBA_BABA total vs D+', paste('num sites >', as.numeric(win_size) / 5)))
abline(h=quantile(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], .99), col='red')

plot(wins$baaa_abaa[wins$num_sites > as.numeric(win_size) / 5], wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], pch=20, col=adjustcolor('gray', .6), xlim=c(0, max(wins$ab[wins$num_sites > as.numeric(win_size) / 5])), xlab='BAAA_ABAA total', ylab='D+', main=c('BAAA_ABAA total vs D+', paste('num sites >', as.numeric(win_size) / 5)))
abline(h=quantile(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], .99), col='red')
par(mfrow=c(1,1))


#traces of D+ windows for each chrom
for (chrom in chroms)
  {
  win_filter <- wins$chrom == chrom & wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff
  plot((wins$start[win_filter] + wins$end[win_filter]) / 2, wins$d_plus[win_filter], type='l', ylim=c(min(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]), max(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff])), xlab='pos', ylab='D+', main=paste('D+', 'window size =', win_size, chrom))
  abline(h=0, col='grey')

  points(gwas_pos$pos[gwas_pos$chrom == chrom], rep(max(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]), sum(gwas_pos$chrom == chrom)), pch=20, col='red') 
  text(min(wins$start[win_filter]), .95 * max(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]), 'GWAS', col='red')
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('red', .2), border=NA)
  }

#scatterplots of D+ vs site pattern differences
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], (wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], pch=20, xlab='D+', ylab='ABBA - BABA', main=paste('D+ vs ABBA - BABA for window size =', win_size))
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], (wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], pch=20, xlab='D+', ylab='BAAA - ABAA', main=paste('D+ vs BAAA - ABAA for window size =', win_size))


#abba - baba vs baaa - abaa
#colored by D+ value
d_plus_col_pal <- adjustcolor(matlab.like(101), .4)
d_plus_cols <- d_plus_col_pal[1 + round(100 * (wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff] - min(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff])) / (max(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]) - min(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff])))]
plot((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], (wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], pch=20, col=d_plus_cols, xlab='ABBA - BABA', ylab='BAAA - ABAA', main=paste('ABBA - BABA vs BAAA - ABAA for window size =', win_size))

abline(v=0, col='gray')
abline(h=0, col='gray')


#plot custom legend
x_len <- max((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]) - min((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]) 
y_len <- max((wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]) - min((wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]) 

x_min <- min((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff])
x_max <- x_min + 0.05 * x_len
y_min <- min((wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]) + 0.7 * y_len
y_max <- y_min + 0.25 * y_len
y_step <- (y_max - y_min) / length(d_plus_col_pal)

for (i in 1:length(d_plus_cols))
  {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=d_plus_col_pal[i], border=d_plus_col_pal[i])}
  
legend_d_plus <- round(seq(min(round(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], 2)), max(round(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], 2)), length.out=5),2)
text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_d_plus)), legend_d_plus, cex=.6)
text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'D+')


#colored by D+ > or < 0
d_plus_binary_cols <- rep(adjustcolor('black', .4), sum(wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff))
d_plus_binary_cols[wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff] < 0] <- adjustcolor('red', .4)
plot((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], (wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], pch=20, col=d_plus_binary_cols, xlab='ABBA - BABA', ylab='BAAA - ABAA', main=paste('ABBA - BABA vs BAAA - ABAA for window size =', win_size))
abline(v=0, col='gray')
abline(h=0, col='gray')
legend('topleft', c('D+ >= 0', 'D+ < 0'), fill=c('black', 'red'), border=NA)


#gwas p vs D+ colored by chrom
chrom_col_pal <- c('red', 'orange', 'green', 'blue', 'purple')
chrom_cols <- c()
for (i in 1:nrow(gwas_pos))
  {chrom_cols <-c(chrom_cols, chrom_col_pal[which(chroms == gwas_pos$chrom[i])])}

plot(log(gwas_pos$p, 10), gwas_pos$d_plus, pch=20, cex=2, col=chrom_cols, xlab='log_10(p)', ylab='D+', main=paste('GWAS p vs D+ for window size =', win_size))
legend('bottomleft', chroms, fill=chrom_col_pal, border=NA)


#D+ vs % abba and baa sites
per_abba <- 100 * (wins$abba / (wins$abba + wins$baba + wins$baaa + wins$abaa))[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]
per_baaa <- 100 * (wins$baaa / (wins$abba + wins$baba + wins$baaa + wins$abaa))[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff]
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], per_abba[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], pch=20, xlab='D+', ylab='% ABBA sites', main=paste('D+ vs % ABBA sites for window size =', win_size))
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], per_baaa[wins$num_sites > as.numeric(win_size) / 5 & wins$ab > ab_cutoff], pch=20, xlab='D+', ylab='% BAAA sites', main=paste('D+ vs % BAAA sites for window size =', win_size))
dev.off()

