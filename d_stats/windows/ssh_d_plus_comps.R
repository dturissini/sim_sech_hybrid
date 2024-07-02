library("RSQLite")  


args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]




#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('ssh_d_plus_comps_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)

d_stats_table <- paste('d_stat_win_', win_size, sep='')
d_stats_no_denis_table <- paste('d_stat_win_', win_size, '_no_denis', sep='')
d_stats_no_denis_unk_table <- paste('d_stat_win_', win_size, '_no_denis_unk', sep='')
d_stats_only_denis_table <- paste('d_stat_win_', win_size, '_only_denis', sep='')

num_sites_cutoff <- as.numeric(win_size) / 50

d_wins <- dbGetQuery(conn, paste("select a.chrom, a.start, a.end, 
                                  a.d_plus, nd.d_plus d_plus_no_denis, ndu.d_plus d_plus_no_denis_unk, od.d_plus d_plus_only_denis,
                                  a.abba, nd.abba abba_no_denis, ndu.abba abba_no_denis_unk, od.abba abba_only_denis,
                                  a.baba, nd.baba baba_no_denis, ndu.baba baba_no_denis_unk, od.baba baba_only_denis,
                                  a.baaa, nd.baaa baaa_no_denis, ndu.baaa baaa_no_denis_unk, od.baaa baaa_only_denis,
                                  a.abaa, nd.abaa abaa_no_denis, ndu.abaa abaa_no_denis_unk, od.abaa abaa_only_denis                                  
                                  from ", d_stats_table, " a, ", d_stats_no_denis_table, " nd, ", d_stats_no_denis_unk_table, " ndu, ", d_stats_only_denis_table, " od
                                  where a.chrom = nd.chrom
                                  and a.chrom = ndu.chrom
                                  and a.chrom = od.chrom
                                  and a.start = nd.start
                                  and a.start = ndu.start
                                  and a.start = od.start
                                  and a.d_plus is not null
                                  and nd.d_plus is not null
                                  and ndu.d_plus is not null
                                  and od.d_plus is not null
                                  and a.num_sites > ", num_sites_cutoff, sep=''))




pdf(pdf_file, height=12, width=12)
par(mfrow=c(2,2))
#D+
plot(d_wins$d_plus, d_wins$d_plus_no_denis, pch=20, col=adjustcolor('gray', .8), xlab='D+ all', ylab='D+ no sech Denis', main='D+ all vs no sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$d_plus, d_wins$d_plus_no_denis_unk, pch=20, col=adjustcolor('gray', .8), xlab='D+ all', ylab='D+ no sech Denis or Unknown', main='D+ all vs no sech Denis or Unknown')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$d_plus, d_wins$d_plus_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='D+ all', ylab='D+ only sech Denis', main='D+ all vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$d_plus_no_denis, d_wins$d_plus_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='D+ no sech Denis', ylab='D+ only sech Denis', main='D+ no sech Denis vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))


#ABBA
plot(d_wins$abba, d_wins$abba_no_denis, pch=20, col=adjustcolor('gray', .8), xlab='ABBA all', ylab='ABBA no sech Denis', main='ABBA all vs no sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$abba, d_wins$abba_no_denis_unk, pch=20, col=adjustcolor('gray', .8), xlab='ABBA all', ylab='ABBA no sech Denis or Unknown', main='ABBA all vs no sech Denis or Unknown')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$abba, d_wins$abba_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='ABBA all', ylab='ABBA only sech Denis', main='ABBA all vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$abba_no_denis, d_wins$abba_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='ABBA no sech Denis', ylab='ABBA only sech Denis', main='ABBA no sech Denis vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))


#BABA
plot(d_wins$baba, d_wins$baba_no_denis, pch=20, col=adjustcolor('gray', .8), xlab='BABA all', ylab='BABA no sech Denis', main='BABA all vs no sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$baba, d_wins$baba_no_denis_unk, pch=20, col=adjustcolor('gray', .8), xlab='BABA all', ylab='BABA no sech Denis or Unknown', main='BABA all vs no sech Denis or Unknown')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$baba, d_wins$baba_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='BABA all', ylab='BABA only sech Denis', main='BABA all vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$baba_no_denis, d_wins$baba_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='BABA no sech Denis', ylab='BABA only sech Denis', main='BABA no sech Denis vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))


#BAAA
plot(d_wins$baaa, d_wins$baaa_no_denis, pch=20, col=adjustcolor('gray', .8), xlab='BAAA all', ylab='BAAA no sech Denis', main='BAAA all vs no sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$baaa, d_wins$baaa_no_denis_unk, pch=20, col=adjustcolor('gray', .8), xlab='BAAA all', ylab='BAAA no sech Denis or Unknown', main='BAAA all vs no sech Denis or Unknown')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$baaa, d_wins$baaa_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='BAAA all', ylab='BAAA only sech Denis', main='BAAA all vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$baaa_no_denis, d_wins$baaa_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='BAAA no sech Denis', ylab='BAAA only sech Denis', main='BAAA no sech Denis vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))


#ABAA
plot(d_wins$abaa, d_wins$abaa_no_denis, pch=20, col=adjustcolor('gray', .8), xlab='ABAA all', ylab='ABAA no sech Denis', main='ABAA all vs no sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$abaa, d_wins$abaa_no_denis_unk, pch=20, col=adjustcolor('gray', .8), xlab='ABAA all', ylab='ABAA no sech Denis or Unknown', main='ABAA all vs no sech Denis or Unknown')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$abaa, d_wins$abaa_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='ABAA all', ylab='ABAA only sech Denis', main='ABAA all vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))

plot(d_wins$abaa_no_denis, d_wins$abaa_only_denis, pch=20, col=adjustcolor('gray', .8), xlab='ABAA no sech Denis', ylab='ABAA only sech Denis', main='ABAA no sech Denis vs only sech Denis')
abline(0, 1, col=adjustcolor('red', .4))
par(mfrow=c(1,1))
dev.off()

