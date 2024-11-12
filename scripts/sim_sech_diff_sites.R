library(RSQLite)
library(colorRamps)
setwd('/proj/matutelb/projects/gwas/genotype_datasets/sech_oa/')


gwas_results_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_gwas.oa_resistance.glm.linear'
diff_sites_db <- "sim_sech_diff_sites.db"

myCon <- dbConnect(dbDriver("SQLite"), diff_sites_db)


diff_sites <- dbGetQuery(myCon, "select chrom, pos, num_sim, num_sech, num_sech_denis, sim_alt_freq, sech_alt_freq, sech_denis_alt_freq
                                 from sim_sech_diff_sites
                                 where sim_alt_freq = 0")


hybrid_site_counts <- dbGetQuery(myCon, "select chrom, pos, count(*) num_lines
                                        from hybrid_sech_alleles
                                        group by chrom, pos")

hybrid_sech_sites <- dbGetQuery(myCon, "select h.chrom, h.pos, hybrid_line, num_sech_alleles
                                        from hybrid_sech_alleles h, sim_sech_diff_sites s
                                        where h.chrom = s.chrom
                                        and h.pos = s.pos
                                        and num_sech >= 20
                                        and num_sim >= 60")
                                        
hybrid_win <- dbGetQuery(myCon, "select win_size, chrom, start, end, hybrid_line, num_sech_sites, total_diff_sites
                                 from hybrid_sech_win")

win_sizes <- sort(unique(hybrid_win$win_size))
hybrid_lines <- sort(unique(hybrid_sech_sites$hybrid_line))                                        
chroms <- unique(diff_sites$chrom)

gwas <- read.table(gwas_results_file, comment.char='', header=T)
names(gwas)[1] <- 'CHROM'

gwas_filtered <- subset(gwas, gwas$P < 1e-7)


per_win_plots <- function(win_size, hybrid_win, diff_sites, gwas_filtered)
  {
  hist(100 * hybrid_win$num_sech_sites[hybrid_win$win_size == win_size] / hybrid_win$total_diff_sites[hybrid_win$win_size == win_size], breaks = 0:101, col='black', xlab='percent sech alleles', ylab='windows', main=c("Percent sech allles in window of total fixed sites", win_size))
  
  par(mar=c(5.1, 7.1, 4.1, 2.1))
  num_cols <- 20
  win_cols <- adjustcolor(matlab.like(num_cols))
  for (i in 1:length(win_cols))
    {win_cols[i] <- adjustcolor(win_cols[i], i / length(win_cols))}
    
  for (chrom in chroms)
    {
    plot(1, type='n', xlim=c(0, max(diff_sites$pos[diff_sites$chrom == chrom])), ylim=c(0.5, length(hybrid_lines) + 0.5), yaxt='n', xlab='pos', ylab='', main=c(paste('Percent hybrid sech alleles; win_size =', win_size), chrom))
    axis(2, at=seq(1:(length(hybrid_lines) + 1)), labels=c(hybrid_lines, 'Top GWAS'), las=1, cex.axis=0.5)
    
    #sech site windows
    for (i in 1:length(hybrid_lines))
      {
      for (j in which(hybrid_win$chrom == chrom & hybrid_win$hybrid_line == hybrid_lines[i] & hybrid_win$win_size == win_size))
        {
        col_i <- min(ceiling(hybrid_win$num_sech_sites[j] / hybrid_win$total_diff_sites[j] * num_cols), num_cols)
        if (col_i > 0)
          {rect(hybrid_win$start[j], i - 0.25, hybrid_win$end[j], i + 0.25, col=win_cols[col_i], border=NA)}
        }
      }


    #draw legend
    rect_start <- -.02 * max(diff_sites$pos[diff_sites$chrom == chrom])
    for (col_k in 1:num_cols)
      {
      rect(rect_start, col_k / 2 - 0.25, 0, col_k / 2  + 0.25, col=win_cols[col_k], border=NA)
      if (col_k %% 5 == 0)
        {text(1.35 * rect_start, col_k / 2 + .25, col_k / num_cols * 100, cex=.5)}
      }
    
    #gwas results
    chrom_i <- which(chroms == chrom)
    points(gwas_filtered$POS[gwas_filtered$CHROM == chrom_i], rep(length(hybrid_lines) + 1, sum(gwas_filtered$CHROM == chrom_i)), pch=20) 
    abline(h=length(hybrid_lines) + 0.5)   
    }
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  }


pdf('sim_sech_diff_sites.pdf', height=8, width=10.5)
hist(diff_sites$num_sim, breaks = seq(0, max(diff_sites$num_sim) + 1), col='black', xlab='simulans lines', ylab='sites', main='Number of simulans lines')
abline(v=60, col='red')

hist(diff_sites$num_sech, breaks = seq(0, max(diff_sites$num_sech) + 1), col='black', xlab='sechellia lines', ylab='sites', main='Number of sechellia lines')
abline(v=20, col='red')

hist(diff_sites$num_sech_denis, breaks = seq(0, max(diff_sites$num_sech_denis) + 1), col='black', xlab='sechellia lines (Denis only)', ylab='sites', main='Number of sechellia lines (Denis only)')
abline(v=5, col='red')


hist(diff_sites$sech_alt_freq[diff_sites$sech_denis_alt_freq == 1], breaks = seq(-.01, 1, .01), col='black', xlab='sechellia freq', ylab='sites', main='Overall sechellia freq where Denis only is fixed')
hist(diff_sites$sech_alt_freq[diff_sites$sech_denis_alt_freq == 1 & diff_sites$sech_alt_freq < 1], breaks = seq(-.01, 1, .01), col='black', xlab='sechellia freq', ylab='sites', main=c('Overall sechellia freq where Denis only is fixed', 'freq < 1'))
hist(diff_sites$num_sech_denis[diff_sites$sech_denis_alt_freq == 1 & diff_sites$sech_alt_freq < 1], breaks = seq(0, max(diff_sites$num_sech_denis) + 1), col='black', xlab='sechellia lines (Denis only)', ylab='sites', main='Number of sechellia lines (Denis only) where Denis is fixed but sechellia poly')

par(mfrow=c(2,1))
for (chrom in chroms)
  {
  hist(diff_sites$pos[diff_sites$chrom == chrom & diff_sites$sech_alt_freq == 1], breaks=seq(0, max(diff_sites$pos[diff_sites$chrom == chrom]) + 10000, 10000), col='black', xlab='pos', ylab='sites', main=c('Fixed sites: sim-sech', chrom))
  hist(diff_sites$pos[diff_sites$chrom == chrom & diff_sites$sech_denis_alt_freq == 1 & diff_sites$num_sech_denis >= 5], breaks=seq(0, max(diff_sites$pos[diff_sites$chrom == chrom]) + 10000, 10000), col='black', xlab='pos', ylab='sites', main=c('Fixed sites: sim-sech (Denis only; >=5 lines)', chrom))
  }
par(mfrow=c(1,1))


hist(hybrid_site_counts$num_lines, breaks=seq(1, max(hybrid_site_counts$num_lines) + 1), col='black', xlab="Hybrid lines", ylab="Number of sites", main="Hybrid lines per sech allele site")


#percent sech sites per window
for (win_size in win_sizes)
  {
  per_win_plots(win_size, hybrid_win, diff_sites, gwas_filtered)
  }
 

#show all sech sites
par(mar=c(5.1, 7.1, 4.1, 2.1))
for (chrom in chroms)
  {
  plot(1, type='n', xlim=c(0, max(diff_sites$pos[diff_sites$chrom == chrom])), ylim=c(0.5, length(hybrid_lines) + 0.5), yaxt='n', xlab='pos', ylab='', main=c('hybrid sech alleles', chrom))
  axis(2, at=seq(1:length(hybrid_lines)), labels=hybrid_lines, las=1, cex.axis=0.5)
  for (i in 1:length(hybrid_lines))
    {
    points(hybrid_sech_sites$pos[hybrid_sech_sites$chrom == chrom & hybrid_sech_sites$hybrid_line == hybrid_lines[i]], rep(i, sum(hybrid_sech_sites$chrom == chrom & hybrid_sech_sites$hybrid_line == hybrid_lines[i])), pch=20, col=adjustcolor('grey', .6))
    }
  }
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


pdf('sim_sech_diff_win_100000.pdf', height=8, width=10.5)
win_size <- 100000
per_win_plots(win_size, hybrid_win, diff_sites, gwas_filtered)
dev.off()

