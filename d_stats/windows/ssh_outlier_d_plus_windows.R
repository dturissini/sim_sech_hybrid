library("RSQLite")  
library(fields)
library(colorRamps)

args <- commandArgs(trailingOnly = TRUE)
win_size <- as.character(args[1])
outlier_type <- args[2]

#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
gwas_snp_db_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db'
anno_db_file <-  paste('/proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sech_oa_anno.db', sep='')
pdf_file <- paste('ssh_outlier_', outlier_type, '_win_', win_size, '.pdf', sep='')
orthodb_file <-  '/proj/matutelb/projects/gwas/orthodb/odb11v0_subsetted_sim_mel.db'
mel_species_id <- '7227_0'

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, paste("attach database '", gwas_snp_db_file, "' as g", sep=''))
dbSendQuery(conn, paste("attach database '", anno_db_file, "' as a", sep=''))
dbSendQuery(conn, paste("attach database '", orthodb_file, "' as o", sep=''))

#load window and gwas data from db
d_stats_table <- paste('d_stat_win_', win_size, sep='')
win_site_table <- paste('outlier_', outlier_type, '_win_sites_', win_size, sep='')
tmp_dsw_table <- paste(outlier_type, '_tmp_dsw_', win_size, sep='')

d_wins <- dbGetQuery(conn, paste("select *
                                  from ", d_stats_table, "
                                  where d_plus is not null", sep=''))

win_sites <- dbGetQuery(conn, paste("select *, 
                                     der_alleles_sim / total_alleles_sim der_freq_sim,
                                     der_alleles_ssh / total_alleles_ssh der_freq_ssh,
                                     der_alleles_sech / total_alleles_sech der_freq_sech 
                                     from ", win_site_table, "
                                     order by chrom, pos", sep=''))


gwas_pos <- dbGetQuery(conn, paste("select dsw_id, g.chrom, pos, p, start, end, d_plus
                                    from gwas_snp g, ", d_stats_table, " w
                                    where g.chrom = w.chrom
                                    and pos between start and end
                                    and p < 1e-7", sep=''))

dbSendQuery(conn, paste("drop table if exists", tmp_dsw_table))
dbSendQuery(conn, paste("create table ", tmp_dsw_table, " as
                         select distinct s.dsw_id, w.chrom, start, end
                         from ", win_site_table, " s, ", d_stats_table, " w
                         where s.dsw_id = w.dsw_id", sep=''))

dbSendQuery(conn, paste("create index idx_se_", outlier_type, win_size, " on ", tmp_dsw_table, "(start, end)", sep=''))
dbSendQuery(conn, paste("create index idx_e_", outlier_type, win_size, " on ", tmp_dsw_table, "(end)", sep=''))

genes <- dbGetQuery(conn, paste("select w.dsw_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from ", tmp_dsw_table, " w, gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and gi.start between w.start and w.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.dsw_id, gi.start, gi.end
                                 union
                                 select w.dsw_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from ", tmp_dsw_table, " w, gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and gi.end between w.start and w.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.dsw_id, gi.start, gi.end
                                 union
                                 select w.dsw_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from gene_info gi, ", tmp_dsw_table, " w, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and w.start between gi.start and gi.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.dsw_id, gi.start, gi.end
                                 union
                                 select w.dsw_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from gene_info gi, ", tmp_dsw_table, " w, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and w.end between gi.start and gi.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.dsw_id, gi.start, gi.end
                                 order by w.dsw_id, gi.start", sep='')) 

dbSendQuery(conn, paste("drop table if exists", tmp_dsw_table))


pdf(pdf_file, height=8, width=10.5)
layout(matrix(c(1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 5, 6, 7, 7, 7, 8, 9, 9, 9, 10), nrow=5, byrow=T))
for (dsw_id in sort(unique(win_sites$dsw_id)))
  {
  #sim
  par(mar=c(0, 4.1, 4.1, 0))
  plot(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sim[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sim[win_sites$dsw_id == dsw_id], pch=20, col='blue', xaxt='n', xlab='', ylim=c(0,1), ylab='Derived freq', main='')
  abline(v=gwas_pos$pos[gwas_pos$dsw_id == dsw_id], col='black')
  if (sum(genes$dsw_id == dsw_id) > 0)
    {
    rect(genes$start[genes$dsw_id == dsw_id], -1, genes$end[genes$dsw_id == dsw_id], 2, col=adjustcolor('gray', .2), border=NA)
    mtext(substr(genes$synonyms[genes$dsw_id == dsw_id], 1, 30), side = 3, line = 0:2, at = (sapply(genes$start[genes$dsw_id == dsw_id], max, min(win_sites$pos[win_sites$dsw_id == dsw_id])) + sapply(genes$end[genes$dsw_id == dsw_id], min, max(win_sites$pos[win_sites$dsw_id == dsw_id]))) / 2, cex=.5)
    }

  par(mar=c(0, 0, 4.1, .1))
  sim_hist <- hist(win_sites$der_alleles_sim[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sim > 0] / win_sites$total_alleles_sim[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sim > 0], breaks=seq(0, 1.01, .01), plot=F)
  plot(1, type='n', xlim=c(0, max(sim_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', main='')
  polygon(c(0, 0, sim_hist$counts, 1), c(0, sim_hist$mids[1], sim_hist$mids, sim_hist$mids[length(sim_hist$mids)]), col='blue', border=NA)  
  legend('topright', 'sim', fill='blue', border=NA, cex=3, bty='n')
  mtext(paste('D+ = ', round(d_wins$d_plus[d_wins$dsw_id == dsw_id], 2), ' (', round(sum(d_wins$d_plus[d_wins$dsw_id == dsw_id] > d_wins$d_plus) / length(d_wins$d_plus), 3), ' percentile)', sep=''), 0, side=3, line=1, at=0, adj=0)
    
  #ssh
  par(mar=c(0, 4.1, 0, 0))
  plot(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_ssh[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_ssh[win_sites$dsw_id == dsw_id], pch=20, col='gold', xaxt='n', xlab='', ylim=c(0,1), ylab='Derived freq', main='')

#  abba_filter <- win_sites$der_freq_sim[win_sites$dsw_id == dsw_id] == 0 & win_sites$der_freq_ssh[win_sites$dsw_id == dsw_id] == 1 & win_sites$der_freq_sech[win_sites$dsw_id == dsw_id] == 1
#  baba_filter <- win_sites$der_freq_sim[win_sites$dsw_id == dsw_id] == 1 & win_sites$der_freq_ssh[win_sites$dsw_id == dsw_id] == 0 & win_sites$der_freq_sech[win_sites$dsw_id == dsw_id] == 1
#  baaa_filter <- win_sites$der_freq_sim[win_sites$dsw_id == dsw_id] == 1 & win_sites$der_freq_ssh[win_sites$dsw_id == dsw_id] == 0 & win_sites$der_freq_sech[win_sites$dsw_id == dsw_id] == 0
#  abaa_filter <- win_sites$der_freq_sim[win_sites$dsw_id == dsw_id] == 0 & win_sites$der_freq_ssh[win_sites$dsw_id == dsw_id] == 1 & win_sites$der_freq_sech[win_sites$dsw_id == dsw_id] == 0
#  
#  points(win_sites$pos[abba_filter], win_sites$der_alleles_ssh[abba_filter], pch='x')
#  points(win_sites$pos[baba_filter], win_sites$der_alleles_ssh[baba_filter], pch='x')
#  points(win_sites$pos[baaa_filter], win_sites$der_alleles_ssh[baaa_filter], pch='O')
#  points(win_sites$pos[abaa_filter], win_sites$der_alleles_ssh[abaa_filter], pch='O')

  abline(v=gwas_pos$pos[gwas_pos$dsw_id == dsw_id], col='black')
  if (sum(genes$dsw_id == dsw_id) > 0)
    {
    rect(genes$start[genes$dsw_id == dsw_id], -1, genes$end[genes$dsw_id == dsw_id], 2, col=adjustcolor('gray', .2), border=NA)
    }

  par(mar=c(0, 0, 0, .1))
  ssh_hist <- hist(win_sites$der_alleles_ssh[win_sites$dsw_id == dsw_id & win_sites$der_alleles_ssh > 0] / win_sites$total_alleles_ssh[win_sites$dsw_id == dsw_id & win_sites$der_alleles_ssh > 0], breaks=seq(0, 1.01, .01), plot=F)
  plot(1, type='n', xlim=c(0, max(ssh_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', main='')
  polygon(c(0, 0, ssh_hist$counts, 1), c(0, ssh_hist$mids[1], ssh_hist$mids, ssh_hist$mids[length(ssh_hist$mids)]), col='gold', border=NA)  
  legend('topright', 'ssh', fill='gold', border=NA, cex=3, bty='n')


  #sech
  par(mar=c(0, 4.1, 0, 0))
  plot(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech[win_sites$dsw_id == dsw_id], pch=20, col='red', xaxt='n', ylim=c(0,1), ylab='Derived freq', main='')
  abline(v=gwas_pos$pos[gwas_pos$dsw_id == dsw_id], col='black')
  if (sum(genes$dsw_id == dsw_id) > 0)
    {
    rect(genes$start[genes$dsw_id == dsw_id], -1, genes$end[genes$dsw_id == dsw_id], 2, col=adjustcolor('gray', .2), border=NA)
    }

  par(mar=c(0, 0, 0, .1))
  sech_hist <- hist(win_sites$der_alleles_sech[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech > 0] / win_sites$total_alleles_sech[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech > 0], breaks=seq(0, 1.01, .01), plot=F)
  plot(1, type='n', xlim=c(0, max(sech_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', main='')
  polygon(c(0, 0, sech_hist$counts, 1), c(0, sech_hist$mids[1], sech_hist$mids, sech_hist$mids[length(sech_hist$mids)]), col='red', border=NA)  
  legend('topright', 'sech', fill='red', border=NA, cex=3, bty='n')


  #sech locations
  par(mar=c(0, 4.1, 0, 0))
  plot(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech_unknown[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech_unknown[win_sites$dsw_id == dsw_id], pch=20, col='gray', xaxt='n', ylim=c(0,1), ylab='Derived freq', main='')

  points(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech_anro[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech_anro[win_sites$dsw_id == dsw_id], pch=20, col='orange')
  points(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech_denis[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech_denis[win_sites$dsw_id == dsw_id], pch=20, col='purple')
  points(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech_ladigue[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech_ladigue[win_sites$dsw_id == dsw_id], pch=20, col='green')
  points(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech_marianne[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech_marianne[win_sites$dsw_id == dsw_id], pch=20, col='darkblue')
  points(win_sites$pos[win_sites$dsw_id == dsw_id], win_sites$der_alleles_sech_praslin[win_sites$dsw_id == dsw_id] / win_sites$total_alleles_sech_praslin[win_sites$dsw_id == dsw_id], pch=20, col='lightblue')


  abline(v=gwas_pos$pos[gwas_pos$dsw_id == dsw_id], col='black')
  if (sum(genes$dsw_id == dsw_id) > 0)
    {
    rect(genes$start[genes$dsw_id == dsw_id], -1, genes$end[genes$dsw_id == dsw_id], 2, col=adjustcolor('gray', .2), border=NA)
    }

  par(mar=c(0, 0, 0, .1))
  sech_unknown_hist <- hist(win_sites$der_alleles_sech_unknown[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_unknown > 0] / win_sites$total_alleles_sech_unknown[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_unknown > 0], breaks=seq(0, 1.01, .01), plot=F)
  sech_anro_hist <- hist(win_sites$der_alleles_sech_anro[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_anro > 0] / win_sites$total_alleles_sech_anro[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_anro > 0], breaks=seq(0, 1.01, .01), plot=F)
  sech_denis_hist <- hist(win_sites$der_alleles_sech_denis[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_denis > 0] / win_sites$total_alleles_sech_denis[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_denis > 0], breaks=seq(0, 1.01, .01), plot=F)
  sech_ladigue_hist <- hist(win_sites$der_alleles_sech_ladigue[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_ladigue > 0] / win_sites$total_alleles_sech_ladigue[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_ladigue > 0], breaks=seq(0, 1.01, .01), plot=F)
  sech_marianne_hist <- hist(win_sites$der_alleles_sech_marianne[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_marianne > 0] / win_sites$total_alleles_sech_marianne[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_marianne > 0], breaks=seq(0, 1.01, .01), plot=F)
  sech_praslin_hist <- hist(win_sites$der_alleles_sech_praslin[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_praslin > 0] / win_sites$total_alleles_sech_praslin[win_sites$dsw_id == dsw_id & win_sites$der_alleles_sech_praslin > 0], breaks=seq(0, 1.01, .01), plot=F)
  plot(1, type='n', xlim=c(0, max(sech_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', xlab='num sites', main='')
  points(sech_unknown_hist$counts, sech_unknown_hist$mids, col='gray', type='l')  
  points(sech_anro_hist$counts, sech_anro_hist$mids, col='orange', type='l')  
  points(sech_denis_hist$counts, sech_denis_hist$mids, col='purple', type='l')  
  points(sech_ladigue_hist$counts, sech_ladigue_hist$mids, col='green', type='l')  
  points(sech_marianne_hist$counts, sech_marianne_hist$mids, col='darkblue', type='l')  
  points(sech_praslin_hist$counts, sech_praslin_hist$mids, col='lightblue', type='l')  
  legend('topright', c('sech - Unk', 'sech - Anro', 'sech - Denis', 'sech - La Digue', 'sech - Marianne', 'sech - Praslin'), fill=c('gray', 'orange', 'purple', 'green', 'darkblue', 'lightblue'), border=NA, cex=1, bty='n')

  



  #ssh freq diffs
  par(mar=c(5.1, 4.1, 0, 0))
  ssh_sim_filter <- win_sites$dsw_id == dsw_id
  ssh_sech_filter <- win_sites$dsw_id == dsw_id
  
  ssh_sim_diff <- abs(win_sites$der_alleles_ssh[ssh_sim_filter] / win_sites$total_alleles_ssh[ssh_sim_filter] - win_sites$der_alleles_sim[ssh_sim_filter] / win_sites$total_alleles_sim[ssh_sim_filter])
  ssh_sech_diff <- abs(win_sites$der_alleles_ssh[ssh_sech_filter] / win_sites$total_alleles_ssh[ssh_sech_filter] - win_sites$der_alleles_sech[ssh_sech_filter] / win_sites$total_alleles_sech[ssh_sech_filter])
  
  plot(win_sites$pos[ssh_sim_filter], ssh_sim_diff, pch=20, col='green', xlab=paste('Pos on', d_wins$chrom[d_wins$dsw_id == dsw_id]), ylim=c(0,1), ylab='Derived freq diff', main='')
  points(win_sites$pos[ssh_sech_filter], ssh_sech_diff, pch=20, col='orange')


  abline(v=gwas_pos$pos[gwas_pos$dsw_id == dsw_id], col='black')
  if (sum(genes$dsw_id == dsw_id) > 0)
    {
    rect(genes$start[genes$dsw_id == dsw_id], -1, genes$end[genes$dsw_id == dsw_id], 2, col=adjustcolor('gray', .2), border=NA)
    }

  par(mar=c(5.1, 0, 0, .1))
  ssh_sim_hist <- hist(ssh_sim_diff, breaks=seq(0, 1.01, .01), plot=F)
  ssh_sech_hist <- hist(ssh_sech_diff, breaks=seq(0, 1.01, .01), plot=F)
  plot(1, type='n', xlim=c(0, max(sech_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', xlab='num sites', main='')
  points(ssh_sim_hist$counts, ssh_sim_hist$mids, col='green', type='l')  
  points(ssh_sech_hist$counts, ssh_sech_hist$mids, col='orange', type='l')  
  legend('topright', c('ssh - sim', 'ssh - sech'), fill=c('green', 'orange'), border=NA, cex=1, bty='n')
  abline(h=mean(ssh_sim_diff), col=adjustcolor('green', .4))
  abline(h=mean(ssh_sech_diff), col=adjustcolor('orange', .4))
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  }  
layout(1)
dev.off()

