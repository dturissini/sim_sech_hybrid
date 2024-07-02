library(RSQLite)
library(colorRamps)
setwd('/proj/matutelb/projects/gwas/genotype_datasets/sech_oa')


freqs_db <- "sim_sech_freq_dist.db"

myCon <- dbConnect(dbDriver("SQLite"), freqs_db)


sim_freqs <- dbGetQuery(myCon, "select sim_alt_freq_bin, sum(num_sites) num_sites
                                from sim_sech_freq_dist
                                where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                group by sim_alt_freq_bin")

sim_freqs_sech_1 <- dbGetQuery(myCon, "select sim_alt_freq_bin, sum(num_sites) num_sites
                                       from sim_sech_freq_dist
                                       where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                       and sech_alt_freq_bin = 1
                                       group by sim_alt_freq_bin")



sim_freqs_ssh_1 <- dbGetQuery(myCon, "select sim_alt_freq_bin, sum(num_sites) num_sites
                                      from sim_sech_freq_dist
                                      where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                      and sim_sech_hybrid_alt_freq_bin = 1
                                      group by sim_alt_freq_bin")

sim_freqs_others_1 <- dbGetQuery(myCon, "select sim_alt_freq_bin, sum(num_sites) num_sites
                                         from sim_sech_freq_dist
                                         where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                         and sech_alt_freq_bin = 1
                                         and sim_sech_hybrid_alt_freq_bin = 1
                                         group by sim_alt_freq_bin")


sech_freqs <- dbGetQuery(myCon, "select sech_alt_freq_bin, sum(num_sites) num_sites
                                 from sim_sech_freq_dist
                                 where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                 group by sech_alt_freq_bin")

ssh_freqs <- dbGetQuery(myCon, "select sim_sech_hybrid_alt_freq_bin, sum(num_sites) num_sites
                                from sim_sech_freq_dist
                                where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                group by sim_sech_hybrid_alt_freq_bin")

ssh_freqs_others_0 <- dbGetQuery(myCon, "select sim_sech_hybrid_alt_freq_bin, sum(num_sites) num_sites
                                from sim_sech_freq_dist
                                where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                and sim_alt_freq_bin = 0 
                                and sech_alt_freq_bin = 0
                                group by sim_sech_hybrid_alt_freq_bin")

ssh_freqs_others_diff <- dbGetQuery(myCon, "select sim_sech_hybrid_alt_freq_bin, sum(num_sites) num_sites
                                            from sim_sech_freq_dist
                                            where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                            and sim_alt_freq_bin = 0 
                                            and sech_alt_freq_bin = 1
                                            group by sim_sech_hybrid_alt_freq_bin")

ssh_dists <- dbGetQuery(myCon, "select sqrt(power(sim_alt_freq_bin - sim_sech_hybrid_alt_freq_bin, 2) + power(sech_alt_freq_bin - sim_sech_hybrid_alt_freq_bin, 2)) ssh_dist, sum(num_sites) num_sites
                                from sim_sech_freq_dist
                                where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                group by sqrt(power(sim_alt_freq_bin - sim_sech_hybrid_alt_freq_bin, 2) + power(sech_alt_freq_bin - sim_sech_hybrid_alt_freq_bin, 2))")


sim_ssh_freqs <- dbGetQuery(myCon, "select sim_alt_freq_bin, sim_sech_hybrid_alt_freq_bin, sum(num_sites) num_sites
                                    from sim_sech_freq_dist
                                    where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                    group by sim_alt_freq_bin, sim_sech_hybrid_alt_freq_bin")

sech_sim_freqs <- dbGetQuery(myCon, "select sech_alt_freq_bin, sim_alt_freq_bin, sum(num_sites) num_sites
                                    from sim_sech_freq_dist
                                    where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                    group by sech_alt_freq_bin, sim_alt_freq_bin")

sech_ssh_freqs <- dbGetQuery(myCon, "select sech_alt_freq_bin, sim_sech_hybrid_alt_freq_bin, sum(num_sites) num_sites
                                    from sim_sech_freq_dist
                                    where not (sim_alt_freq_bin = 0 and sech_alt_freq_bin = 0 and sim_sech_hybrid_alt_freq_bin = 0)
                                    group by sech_alt_freq_bin, sim_sech_hybrid_alt_freq_bin")




chroms <- data.frame(cbind(c('2L', '2R', '3L', '3R', '4', 'X'), c(23857595, 22319025, 23399903, 28149585, 1146867, 22032822)))
names(chroms) <- c('chrom', 'len')
chroms$len <- as.numeric(chroms$len) 


bin_size <- .05
pdf('sim_sech_freqs_all.pdf', height=8, width=10.5)
#allele freq hists
par(mfrow=c(3,1))
plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(sim_freqs$num_sites)), xlab='alt allele freq', ylab='sites', main='sim alt allele freq')
for (bin_start in sim_freqs$sim_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, sim_freqs$num_sites[sim_freqs$sim_alt_freq_bin == bin_start], col='black')}
  
plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(sech_freqs$num_sites)), xlab='alt allele freq', ylab='sites', main='sech alt allele freq')
for (bin_start in sech_freqs$sech_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, sech_freqs$num_sites[sech_freqs$sech_alt_freq_bin == bin_start], col='black')}

plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(ssh_freqs$num_sites)), xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq')
for (bin_start in ssh_freqs$sim_sech_hybrid_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, ssh_freqs$num_sites[ssh_freqs$sim_sech_hybrid_alt_freq_bin == bin_start], col='black')}



plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(sim_freqs_sech_1$num_sites)), xlab='alt allele freq', ylab='sites', main=c('sim alt allele freq', 'sech alt freq = 1'))
for (bin_start in sim_freqs_sech_1$sim_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, sim_freqs_sech_1$num_sites[sim_freqs_sech_1$sim_alt_freq_bin == bin_start], col='black')}

plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(sim_freqs_ssh_1$num_sites)), xlab='alt allele freq', ylab='sites', main=c('sim alt allele freq', 'sim_sech_hybrid alt freq = 1'))
for (bin_start in sim_freqs_ssh_1$sim_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, sim_freqs_ssh_1$num_sites[sim_freqs_ssh_1$sim_alt_freq_bin == bin_start], col='black')}

plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(sim_freqs_others_1$num_sites)), xlab='alt allele freq', ylab='sites', main=c('sim alt allele freq', 'sech and sim_sech_hybrid alt freqs = 1'))
for (bin_start in sim_freqs_others_1$sim_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, sim_freqs_others_1$num_sites[sim_freqs_others_1$sim_alt_freq_bin == bin_start], col='black')}



par(mfrow=c(4,1))
plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(ssh_freqs$num_sites)), xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq')
for (bin_start in ssh_freqs$sim_sech_hybrid_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, ssh_freqs$num_sites[ssh_freqs$sim_sech_hybrid_alt_freq_bin == bin_start], col='black')}

plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(ssh_freqs_others_0$num_sites)), xlab='alt allele freq', ylab='sites', main=c('sim_sech_hybrid alt allele freq', 'sim and sech alt freqs = 0'))
for (bin_start in ssh_freqs_others_0$sim_sech_hybrid_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, ssh_freqs_others_0$num_sites[ssh_freqs_others_0$sim_sech_hybrid_alt_freq_bin == bin_start], col='black')}

plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(ssh_freqs_others_diff$num_sites)), xlab='alt allele freq', ylab='sites', main=c('sim_sech_hybrid alt allele freq', 'alt freqs: sim = 0; sech = 1'))
for (bin_start in ssh_freqs_others_diff$sim_sech_hybrid_alt_freq_bin)
  {rect(bin_start, 0, bin_start + bin_size, ssh_freqs_others_diff$num_sites[ssh_freqs_others_diff$sim_sech_hybrid_alt_freq_bin == bin_start], col='black')}

plot(1, type='n', xlim=c(0, 1.05), ylim=c(0, max(ssh_freqs_others_diff$num_sites[ssh_freqs_others_diff$sim_sech_hybrid_alt_freq_bin != 0])), xlab='alt allele freq', ylab='sites', main=c('sim_sech_hybrid alt allele freq', 'alt freqs: sim = 0; sech = 1'))
for (bin_start in ssh_freqs_others_diff$sim_sech_hybrid_alt_freq_bin[ssh_freqs_others_diff$sim_sech_hybrid_alt_freq_bin != 0])
  {rect(bin_start, 0, bin_start + bin_size, ssh_freqs_others_diff$num_sites[ssh_freqs_others_diff$sim_sech_hybrid_alt_freq_bin == bin_start], col='black')}
par(mfrow=c(1,1))



#define color vectors
freq_cols = c('black', matlab.like(round(1 / bin_size)))
bins <- seq(0, 1, bin_size)
sim_ssh_freq_dist <- matrix(sim_ssh_freqs$num_sites, nrow=21, ncol=21, byrow=T)
sech_sim_freq_dist <- matrix(sech_sim_freqs$num_sites, nrow=21, ncol=21, byrow=T)
sech_ssh_freq_dist <- matrix(sech_ssh_freqs$num_sites, nrow=21, ncol=21, byrow=T)

image(bins, bins, sim_ssh_freq_dist, col=freq_cols, xlab='sim freq', ylab='sim_sech_hybrid freq', main='sim vs sim_sech_hybrid alt freq: all sites')
image(bins, bins, log(sim_ssh_freq_dist, 10), col=freq_cols, xlab='sim freq', ylab='sim_sech_hybrid freq', main=c('sim vs sim_sech_hybrid alt freq: all sites', 'log-scaled'))

image(bins, bins, sech_sim_freq_dist, col=freq_cols, xlab='sech freq', ylab='sim freq', main='sech vs sim alt freq: all sites')
image(bins, bins, log(sech_sim_freq_dist, 10), col=freq_cols, xlab='sech freq', ylab='sim freq', main=c('sech vs sim alt freq: all sites', 'log-scaled'))

image(bins, bins, sech_ssh_freq_dist, col=freq_cols, xlab='sech freq', ylab='sim_sech_hybrid freq', main='sech vs sim_sech_hybrid alt freq: all sites')
image(bins, bins, log(sech_ssh_freq_dist, 10), col=freq_cols, xlab='sech freq', ylab='sim_sech_hybrid freq', main=c('sech vs sim_sech_hybrid alt freq: all sites', 'log-scaled'))
dev.off()




#color by ssh distance
dist_cols <- c('black', matlab.like(100))


########     wortk out hist based on counts for range of values

plot(ssh_dists$ssh_dist, ssh_dists$num_sites, type='l', xlim=c(0, max(ssh_dists$ssh_dist)), ylim=c(0, max(ssh_dists$num_sites)), xlab='sim_sech_hybrid freq distance', ylab='sites', main=c('sim_sech_hybrid alt allele freq distance', 'distance = sqrt((sim - hybrid)^2 + (sech - hybrid)^2)'))
##for (i in 1:length(ssh_hist$counts))
##  {rect(ssh_hist$breaks[i], 0, ssh_hist$breaks[i + 1], ssh_hist$counts[i], col=dist_cols[1 + round(100 * (ssh_hist$breaks[i] + ssh_hist$breaks[i + 1]) / 2 / max(ssh_dist))], border=NA)}



#####par(mfrow=c(2,1))
#####hist(freqs$sim_sech_hybrid_alt_freq, breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main='sim_sech_hybrid alt allele freq')
#####hist(freqs$sim_sech_hybrid_alt_freq[freqs$sim_alt_freq == 1 & freqs$sech_alt_freq == 1], breaks=seq(0, 1.01, .01), col='black', xlab='alt allele freq', ylab='sites', main=c('sim_sech_hybrid alt allele freq', 'both sim and sech alt allele freqs = 1'))
#####par(mfrow=c(1,1))
dev.off()
