library(colorRamps)
setwd('/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression')

pca_raw <- read.table('sim_sech_outgroup_mel_pca.eigenvec', header=T, comment.char='')
samples <- read.table('d_stats/sim_sech_hybrid_ripta_meta.txt', header=F, sep="\t")
names(samples) <- c('sample_id', 'species')

pca <- merge(pca_raw, samples, by.x='IID', by.y= 'sample_id', all.x=T)

cols <- rep('black', nrow(pca))
cols[which(pca$species == 'sechellia')] <- 'red'
cols[which(pca$species == 'sim_sech_hybrid')] <- 'green'
cols[which(pca$species == 'melanogaster')] <- 'purple'


pdf('ssh_mel_out_pca.pdf', height=8, width=10.5)
plot(pca$PC1, pca$PC2,  pch=20, col=cols, xlab='PC1', ylab='PC2',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC2, pca$PC3,  pch=20, col=cols, xlab='PC2', ylab='PC3',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC3, pca$PC4,  pch=20, col=cols, xlab='PC3', ylab='PC4',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC4, pca$PC5,  pch=20, col=cols, xlab='PC4', ylab='PC5',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC5, pca$PC6,  pch=20, col=cols, xlab='PC5', ylab='PC6',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC6, pca$PC7,  pch=20, col=cols, xlab='PC6', ylab='PC7',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC7, pca$PC8,  pch=20, col=cols, xlab='PC7', ylab='PC8',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC8, pca$PC9,  pch=20, col=cols, xlab='PC8', ylab='PC9',  main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)

plot(pca$PC9, pca$PC10, pch=20, col=cols, xlab='PC9', ylab='PC10', main='')
legend('bottomright', c('sim', 'sech', 'sim_sech_hybrid', 'melanogaster'), fill=c('black', 'red', 'green', 'purple'), border=NA)
dev.off()


