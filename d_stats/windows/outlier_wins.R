library("RSQLite")  

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- as.character(args[1])
outlier_type <- args[2]
pop_str <- args[3]


#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
gwas_snp_db_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db'
anno_db_file <-  paste('/proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sech_oa_anno.db', sep='')
pdf_file <- paste('outlier_win_', outlier_type, '_', win_size, '_', pop_str, '.pdf', sep='')
orthodb_file <-  '/proj/matutelb/projects/gwas/orthodb/odb11v0_subsetted_sim_mel.db'
mel_species_id <- '7227_0'

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, paste("attach database '", gwas_snp_db_file, "' as g", sep=''))
dbSendQuery(conn, paste("attach database '", anno_db_file, "' as a", sep=''))
dbSendQuery(conn, paste("attach database '", orthodb_file, "' as o", sep=''))

#define db tables
d_stats_table <- paste('d_stat_win_', win_size, '_', pop_str, sep='')
win_site_table <- paste('outlier_', outlier_type, '_win_sites_', win_size, '_', pop_str, sep='')
tmp_dsw_table <- paste(outlier_type, '_tmp_dsw_', win_size, sep='')

#load window and gwas data from db
d_wins <- dbGetQuery(conn, paste("select *
                                  from ", d_stats_table, "
                                  where d_plus is not null", sep=''))

win_sites <- dbGetQuery(conn, paste("select *
                                     from ", win_site_table, "
                                     where der_alleles > 0
                                     order by chrom, pos", sep=''))

#get pops and colors for plotting
pop_cols <- dbGetQuery(conn, paste("select pop, short_desc, col
                                     from pop_cols
                                     order by pop", sep=''))

gwas_pos <- dbGetQuery(conn, paste("select win_id, g.chrom, pos, p, start, end, d_plus
                                    from gwas_snp g, ", d_stats_table, " w
                                    where g.chrom = w.chrom
                                    and pos between start and end
                                    and p < 1e-7", sep=''))


#make temporary table of outlier windows to speed up queries identifying overlapping homologous melanogaster genes
dbSendQuery(conn, paste("drop table if exists", tmp_dsw_table))
dbSendQuery(conn, paste("create table ", tmp_dsw_table, " as
                         select distinct s.win_id, w.chrom, start, end
                         from ", win_site_table, " s, ", d_stats_table, " w
                         where s.win_id = w.win_id", sep=''))

dbSendQuery(conn, paste("create index idx_se_", outlier_type, win_size, " on ", tmp_dsw_table, "(start, end)", sep=''))
dbSendQuery(conn, paste("create index idx_e_", outlier_type, win_size, " on ", tmp_dsw_table, "(end)", sep=''))

#get genes overlapping windows and identify mel orthologs
#there are four distinct queries to optimize index use for the four different ways a gene and window can overlap. 
#Some genes will be selected in multiple queries but the use of unions will remove duplicates
#1) gene start within window
#2) gene end within window
#3) window start within gene
#4) window end within gene
genes <- dbGetQuery(conn, paste("select w.win_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from ", tmp_dsw_table, " w, gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and gi.start between w.start and w.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.win_id, gi.start, gi.end
                                 union
                                 select w.win_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from ", tmp_dsw_table, " w, gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and gi.end between w.start and w.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.win_id, gi.start, gi.end
                                 union
                                 select w.win_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from gene_info gi, ", tmp_dsw_table, " w, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and w.start between gi.start and gi.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.win_id, gi.start, gi.end
                                 union
                                 select w.win_id, gi.start, gi.end, group_concat(distinct synonyms order by synonyms) synonyms
                                 from gene_info gi, ", tmp_dsw_table, " w, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                 where w.chrom = gi.chrom
                                 and w.end between gi.start and gi.end
                                 and gi.ncbi_gene_id = external_id
                                 and gx.orthodb_gene_id = og.orthodb_gene_id
                                 and og.orthogroup_id = og2.orthogroup_id
                                 and og.orthodb_gene_id != og2.orthodb_gene_id
                                 and og2.orthodb_gene_id = g.orthodb_gene_id
                                 and orthodb_species_id = '", mel_species_id, "'
                                 group by w.win_id, gi.start, gi.end
                                 order by w.win_id, gi.start", sep='')) 

dbSendQuery(conn, paste("drop table if exists", tmp_dsw_table))

#define pop vectors
species_pops <- c('sim', 'sech', 'ssh')
sech_pops <- c('sechanro', 'sechden', 'sechlad', 'sechmari', 'sechpras', 'sechunk')



#plot derived allele freqs for each outlier window
pdf(pdf_file, height=8, width=10.5)
layout(matrix(c(1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 5, 6, 7, 7, 7, 8), nrow=4, byrow=T))
for (win_id in sort(unique(win_sites$win_id)))
  {
  # derived allele freqs
  par(mar=c(0, 4.1, 4.1, 0))
  counter <- 0
  for (pop in species_pops)
    {  
    counter <- counter + 1
    
    plot(win_sites$pos[win_sites$win_id == win_id & win_sites$pop == pop], win_sites$der_alleles[win_sites$win_id == win_id & win_sites$pop == pop] / win_sites$total_alleles[win_sites$win_id == win_id & win_sites$pop == pop], pch=20, col=pop_cols$col[pop_cols$pop == pop], xaxt='n', xlab='', ylim=c(0,1), ylab='Derived freq', main='')
    abline(v=gwas_pos$pos[gwas_pos$win_id == win_id], col='black')
    if (sum(genes$win_id == win_id) > 0)
      {
      rect(genes$start[genes$win_id == win_id], -1, genes$end[genes$win_id == win_id], 2, col=adjustcolor('gray', .2), border=NA)
      if (counter == 1)
        {mtext(substr(genes$synonyms[genes$win_id == win_id], 1, 30), side = 3, line = 0:2, at = (sapply(genes$start[genes$win_id == win_id], max, min(win_sites$pos[win_sites$win_id == win_id])) + sapply(genes$end[genes$win_id == win_id], min, max(win_sites$pos[win_sites$win_id == win_id]))) / 2, cex=.5)}
      }
    
    if (counter == 1)
      {
      par(mar=c(0, 0, 4.1, .1))
      } else {
      par(mar=c(0, 0, 0, .1))
      }
      
    pop_hist <- hist(win_sites$der_alleles[win_sites$win_id == win_id & win_sites$pop == pop] / win_sites$total_alleles[win_sites$win_id == win_id & win_sites$pop == pop], breaks=seq(0, 1.01, .01), plot=F)
    plot(1, type='n', xlim=c(0, max(pop_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', main='')
    polygon(c(0, 0, pop_hist$counts, 1), c(0, pop_hist$mids[1], pop_hist$mids, pop_hist$mids[length(pop_hist$mids)]), col=pop_cols$col[pop_cols$pop == pop], border=NA)  
    legend('topright', pop, fill=pop_cols$col[pop_cols$pop == pop], border=NA, cex=3, bty='n')
    if (counter == 1)
      {mtext(paste('D+ = ', round(d_wins$d_plus[d_wins$win_id == win_id], 2), ' (', round(sum(d_wins$d_plus[d_wins$win_id == win_id] > d_wins$d_plus) / length(d_wins$d_plus), 3), ' percentile)', sep=''), 0, side=3, line=1, at=0, adj=0)}
      
    par(mar=c(0, 4.1, 0, 0))
    }

  #sech pops
  par(mar=c(5.1, 4.1, 0, 0))
  plot(1, type='n', xlim=range(win_sites$pos[win_sites$win_id == win_id]), ylim=c(0,1), ylab='Derived freq', xlab=paste('Position on ', win_sites$chrom[win_sites$win_id == win_id], sep=''), main='')
  abline(v=gwas_pos$pos[gwas_pos$win_id == win_id], col='black')
  if (sum(genes$win_id == win_id) > 0)
    {rect(genes$start[genes$win_id == win_id], -1, genes$end[genes$win_id == win_id], 2, col=adjustcolor('gray', .2), border=NA)}

  for (pop in sech_pops)
    {    
    points(win_sites$pos[win_sites$win_id == win_id & win_sites$pop == pop], win_sites$der_alleles[win_sites$win_id == win_id & win_sites$pop == pop] / win_sites$total_alleles[win_sites$win_id == win_id & win_sites$pop == pop], pch=20, col=pop_cols$col[pop_cols$pop == pop])        
    }
    
  par(mar=c(5.1, 0, 0, .1))
  plot(1, type='n', xlim=c(0, max(pop_hist$counts)), ylim=c(0, 1), xaxt='n', yaxt='n', xlab='num sites', main='')
  for (pop in sech_pops)
    {    
    pop_hist <- hist(win_sites$der_alleles[win_sites$win_id == win_id & win_sites$pop == pop] / win_sites$total_alleles[win_sites$win_id == win_id & win_sites$pop == pop], breaks=seq(0, 1.01, .01), plot=F)
    points(pop_hist$counts, pop_hist$mids, col=pop_cols$col[pop_cols$pop == pop], type='l')  
    }
  legend('topright', pop_cols$short_desc[pop_cols$pop %in% sech_pops], fill=pop_cols$col[pop_cols$pop %in% sech_pops], border=NA, cex=1, bty='n')
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  }  
layout(1)
dev.off()

