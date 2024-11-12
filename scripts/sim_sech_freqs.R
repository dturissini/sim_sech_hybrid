library(RSQLite)
library(colorRamps)
setwd('/proj/matutelb/projects/gwas/genotype_datasets/sech_oa')


freqs_db <- "sim_sech_freqs.db"

myCon <- dbConnect(dbDriver("SQLite"), freqs_db)


freqs <- dbGetQuery(myCon, "select chrom, pos, num_sim, num_sech, num_sech_denis, num_sim_sech_hybrid, sim_alt_freq, sech_alt_freq, sech_denis_alt_freq, sim_sech_hybrid_alt_freq
                            from sim_sech_freqs 
                            where num_sim > 60 
                            and num_sech > 20 
                            and num_sim_sech_hybrid > 10")


chroms <- data.frame(cbind(c('2L', '2R', '3L', '3R', '4', 'X'), c(23857595, 22319025, 23399903, 28149585, 1146867, 22032822)))
names(chroms) <- c('chrom', 'len')
chroms$len <- as.numeric(chroms$len) 



pdf('sim_sech_freqs.pdf', height=8, width=10.5)
#allele freq hists
par(mfrow=c(3,1))
hist(freqs$sim_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim alt allele freq')
hist(freqs$sech_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sech alt allele freq')
hist(freqs$sim_sech_hybrid_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq')
par(mfrow=c(1,1))


#color by sim_sech_hybrid freq
freq_cols = c('black', matlab.like(100))
#plot(freqs$sim_alt_freq, freqs$sech_alt_freq, pch=20, xlab='sim alt allele freq', ylab='sech alt allele freq', main='sim_sech_hybrid freq coloring', col=freq_cols[1 + round(freqs$sim_sech_hybrid_alt_freq * 100)])
#legend('topleft', legend=seq(0, 1, 0.1), fill=freq_cols[1 + seq(0, 100, 10)], border=NA)


#color by ssh distance
ssh_dist <- sqrt((freqs$sim_alt_freq - freqs$sim_sech_hybrid_alt_freq)^2 + (freqs$sech_alt_freq - freqs$sim_sech_hybrid_alt_freq)^2)
ssh_hist <- hist(ssh_dist, breaks=seq(0, max(ssh_dist) + .05, .05), plot=F)
dist_cols <- c('black', matlab.like(100))

plot(1, type='n', xlim=c(0, max(ssh_hist$breaks)), ylim=c(0, max(ssh_hist$counts)), xlab='sim_sech_hybrid freq distance', ylab='sites', main=c('sim_sech_hybrid alt allele freq distance', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'))
for (i in 1:length(ssh_hist$counts))
  {rect(ssh_hist$breaks[i], 0, ssh_hist$breaks[i + 1], ssh_hist$counts[i], col=dist_cols[1 + round(100 * (ssh_hist$breaks[i] + ssh_hist$breaks[i + 1]) / 2 / max(ssh_dist))], border=NA)}
#
#plot(freqs$sim_alt_freq, freqs$sech_alt_freq, pch=20, xlab='sim alt allele freq', ylab='sech alt allele freq', main=c('sim_sech_hybrid distance coloring', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'), col=dist_cols[1 + round(100 * ssh_dist / max(ssh_dist))])
#legend('topleft', legend=round(seq(0, max(ssh_dist), length.out=10), 2), fill=dist_cols[1 + round(100 * seq(0, max(ssh_dist), length.out=10) / max(ssh_dist))], border=NA)
#

par(mfrow=c(2,1))
hist(freqs$sim_sech_hybrid_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq')
hist(freqs$sim_sech_hybrid_alt_freq[freqs$sim_alt_freq == 0 & freqs$sech_alt_freq == 0], breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main=c('sim_sech_hybrid alt allele freq', 'both sim and sech alt allele freqs = 0'))
par(mfrow=c(1,1))
dev.off()

tiff('sim_sech_hybrid_freq_vs_dist.tiff', units="in", width=8, height=5, res=300)
freq_cols <- adjustcolor(freq_cols, .3)
plot(freqs$sim_sech_hybrid_alt_freq, ssh_dist, pch=20, xlab='sim_sech_hybrid alt allele freq', ylab='sim_sech_hybrid freq distance', xlim=c(0,1), col=freq_cols[1 + round(freqs$sech_alt_freq * 100)], main=c('sim_sech_hybrid allele freq vs distance colored by sech freq', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'))
#legend('bottomleft', legend=seq(0, 1, 0.1), fill=freq_cols[1 + seq(0, 100, 10)], border=NA)
dev.off()
