library(RSQLite)
library(colorRamps)
setwd('/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results')


gwas_freqs_db <- "sech_oa_sim_only_pca_gwas_allele_freqs.db"

myCon <- dbConnect(dbDriver("SQLite"), gwas_freqs_db)


gwas_freqs <- dbGetQuery(myCon, "select chrom, pos, gwas_p, num_sim, num_sech, num_sech_denis, num_sim_sech_hybrid, sim_alt_freq, sech_alt_freq, sech_denis_alt_freq, sim_sech_hybrid_alt_freq
                                 from gwas_sim_sech_freqs")


chroms <- data.frame(cbind(c('2L', '2R', '3L', '3R', '4', 'X'), c(23857595, 22319025, 23399903, 28149585, 1146867, 22032822)))
names(chroms) <- c('chrom', 'len')
chroms$len <- as.numeric(chroms$len) 



pdf('sech_oa_sim_only_pca_gwas_allele_freqs.pdf', height=8, width=10.5)
#allele freq hists
par(mfrow=c(3,1))
hist(gwas_freqs$sim_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim alt allele freq: top gwas sites')
hist(gwas_freqs$sech_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sech alt allele freq: top gwas sites')
hist(gwas_freqs$sim_sech_hybrid_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq: top gwas sites')
par(mfrow=c(1,1))


#color by sim_sech_hybrid freq
freq_cols = c('black', matlab.like(100))
plot(gwas_freqs$sim_alt_freq, gwas_freqs$sech_alt_freq, pch=20, xlab='sim alt allele freq', ylab='sech alt allele freq', main='sim_sech_hybrid freq coloring for top gwas snps', col=freq_cols[1 + round(gwas_freqs$sim_sech_hybrid_alt_freq * 100)])
legend('bottomright', legend=seq(0, 1, 0.1), fill=freq_cols[1 + seq(0, 100, 10)], border=NA)


#color by ssh distance
ssh_dist <- sqrt((gwas_freqs$sim_alt_freq - gwas_freqs$sim_sech_hybrid_alt_freq)^2 + (gwas_freqs$sech_alt_freq - gwas_freqs$sim_sech_hybrid_alt_freq)^2)
ssh_hist <- hist(ssh_dist, breaks=seq(0, max(ssh_dist) + .05, .05), plot=F)
dist_cols <- c('black', matlab.like(100))

plot(1, type='n', xlim=c(0, max(ssh_hist$breaks)), ylim=c(0, max(ssh_hist$counts)), xlab='sim_sech_hybrid freq distance', ylab='sites', main=c('sim_sech_hybrid alt allele freq distance for top gwas sites', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'))
for (i in 1:length(ssh_hist$counts))
  {rect(ssh_hist$breaks[i], 0, ssh_hist$breaks[i + 1], ssh_hist$counts[i], col=dist_cols[1 + round(100 * (ssh_hist$breaks[i] + ssh_hist$breaks[i + 1]) / 2 / max(ssh_dist))], border=NA)}

plot(gwas_freqs$sim_alt_freq, gwas_freqs$sech_alt_freq, pch=20, xlab='sim alt allele freq', ylab='sech alt allele freq', main=c('sim_sech_hybrid distance coloring for top gwas snps', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'), col=dist_cols[1 + round(100 * ssh_dist / max(ssh_dist))])
legend('bottomright', legend=round(seq(0, max(ssh_dist), length.out=10), 2), fill=dist_cols[1 + round(100 * seq(0, max(ssh_dist), length.out=10) / max(ssh_dist))], border=NA)

plot(gwas_freqs$sim_sech_hybrid_alt_freq, ssh_dist, pch=20, xlab='sim_sech_hybrid alt allele freq', ylab='sim_sech_hybrid freq distance', xlim=c(0,1), col=freq_cols[1 + round(gwas_freqs$sech_alt_freq * 100)], main=c('sim_sech_hybrid allele freq vs distance colored by sech freq', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'))
legend('bottomleft', legend=seq(0, 1, 0.1), fill=freq_cols[1 + seq(0, 100, 10)], border=NA)

par(mfrow=c(2,1))
hist(gwas_freqs$sim_sech_hybrid_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq: top gwas sites')
hist(gwas_freqs$sim_sech_hybrid_alt_freq[gwas_freqs$sim_alt_freq == 0 & gwas_freqs$sech_alt_freq == 0], breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main=c('sim_sech_hybrid alt allele freq: top gwas sites', 'both sim and sech alt allele freqs = 0'))
par(mfrow=c(1,1))



#manhattan plot
plot(1, type='n', xlim=c(0, sum(chroms$len)), ylim=c(0, max(-log(gwas_freqs$gwas_p, 10))), xaxt='n', xlab='', ylab='-log10(p)', main='gwas with sim_sech_hybrid distance coloring')
legend('topleft', legend=round(seq(0, max(ssh_dist), length.out=10), 2), fill=dist_cols[1 + round(100 * seq(0, max(ssh_dist), length.out=10) / max(ssh_dist))], border=NA)
chrom_offset <- 0
for (chrom in chroms$chrom)
  {
  if (chrom_offset > 0) 
    {abline(v=chrom_offset, col='darkgrey')}
  
  points(gwas_freqs$pos[gwas_freqs$chrom == chrom] + chrom_offset, -log(gwas_freqs$gwas_p[gwas_freqs$chrom == chrom], 10), pch=20, col=dist_cols[1 + round(100 * ssh_dist / max(ssh_dist))])
  
  text(chrom_offset + chroms$len[chroms$chrom == chrom] / 2, 3, chrom, cex=2)
  chrom_offset <- chrom_offset + chroms$len[chroms$chrom == chrom]
  }


plot(1, type='n', xlim=c(0, sum(chroms$len)), ylim=c(0, max(-log(gwas_freqs$gwas_p, 10))), xaxt='n', xlab='', ylab='-log10(p)', main='gwas with sech freq coloring')
legend('topleft', legend=seq(0, 1, 0.1), fill=freq_cols[1 + seq(0, 100, 10)], border=NA)
chrom_offset <- 0
for (chrom in chroms$chrom)
  {
  if (chrom_offset > 0) 
    {abline(v=chrom_offset, col='darkgrey')}

  points(gwas_freqs$pos[gwas_freqs$chrom == chrom] + chrom_offset, -log(gwas_freqs$gwas_p[gwas_freqs$chrom == chrom], 10), pch=20, col=freq_cols[1 + round(gwas_freqs$sech_alt_freq * 100)])

  text(chrom_offset + chroms$len[chroms$chrom == chrom] / 2, 3, chrom, cex=2)
  chrom_offset <- chrom_offset + chroms$len[chroms$chrom == chrom]
  }
dev.off()


