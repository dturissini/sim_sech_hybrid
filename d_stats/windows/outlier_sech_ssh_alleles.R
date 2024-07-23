library("RSQLite")  

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- as.character(args[1])
pop_str <- as.character(args[2])

#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
gwas_snp_db_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db'
anno_db_file <-  paste('/proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sech_oa_anno.db', sep='')
pdf_file <- paste('outlier_d_plus_sech_ssh_alleles_', win_size, '_', pop_str, '.pdf', sep='')
orthodb_file <-  '/proj/matutelb/projects/gwas/orthodb/odb11v0_subsetted_sim_mel.db'
mel_species_id <- '7227_0'

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, paste("attach database '", gwas_snp_db_file, "' as g", sep=''))
dbSendQuery(conn, paste("attach database '", anno_db_file, "' as a", sep=''))
dbSendQuery(conn, paste("attach database '", orthodb_file, "' as o", sep=''))

#load window and gwas data from db
d_stats_table <- paste('d_stat_win_', win_size, '_', pop_str, sep='')
win_site_table <- paste('outlier_d_plus_win_sites_', win_size, '_', pop_str, sep='')
allele_table_sech <- paste('outlier_d_plus_win_alleles_sech_', win_size, '_', pop_str, sep='')
allele_table_ssh <- paste('outlier_d_plus_win_alleles_ssh_', win_size, '_', pop_str, sep='')
allele_dist_table <- paste('outlier_d_plus_win_sech_adj_allele_dist_', win_size, '_', pop_str, sep='')
tmp_dsw_table <- paste('outlier_tmp_dsw_', win_size, sep='')

#define pop str to use in db queries
pop_sql_str <- "'sechanro', 'sechden', 'sechlad', 'sechmari', 'sechpras', 'sechunk', 'sshlad', 'sshmahe', 'sshanro'"



#get win data from db
d_wins <- dbGetQuery(conn, paste("select *
                                  from ", d_stats_table, "
                                  where d_plus is not null", sep=''))

win_sites <- dbGetQuery(conn, paste("select *
                                     from ", win_site_table, "
                                     where 1.0 * der_alleles / total_alleles > 0
                                     and 1.0* der_alleles / total_alleles < 1 
                                     and pop = 'sech'
                                     order by chrom, pos", sep=''))



                                     #get samples and pops from db 
samples <- dbGetQuery(conn, paste("select pop, sample_id
                                   from sample_pop_link
                                   where pop in (", pop_sql_str, ")
                                   order by pop, sample_id", sep=''))

#get pops and colors for plotting
sech_ssh_pops <- dbGetQuery(conn, paste("select pop, short_desc, col
                                         from pop_cols
                                         where pop in (", pop_sql_str, ")
                                         order by pop", sep=''))


#get gwas data from db
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

dbSendQuery(conn, paste("create index idx_se_outlier", win_size, " on ", tmp_dsw_table, "(start, end)", sep=''))
dbSendQuery(conn, paste("create index idx_e_outlier", win_size, " on ", tmp_dsw_table, "(end)", sep=''))

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

#drop tmp table
dbSendQuery(conn, paste("drop table if exists", tmp_dsw_table))


win_size <- as.numeric(win_size)






pdf(pdf_file, height=8, width=10.5)
layout(matrix(c(1, 2, 2, 2, 2), nrow=5, byrow=T))
for (win_id in sort(unique(win_sites$win_id)))
  {
  cat(win_id, as.character(Sys.time()), "\n")
  #ssh freqs
  par(mar=c(0, 10.1, 4.1, 1.1))
  plot(win_sites$pos[win_sites$win_id == win_id], win_sites$der_alleles[win_sites$win_id == win_id] / win_sites$total_alleles[win_sites$win_id == win_id], pch=20, col='red', xaxt='n', ylim=c(0,1), ylab='Der sech freq', main='')
  abline(v=gwas_pos$pos[gwas_pos$win_id == win_id], col='black')

  mtext(paste('D+ = ', round(d_wins$d_plus[d_wins$win_id == win_id], 2), sep=''), 0, cex=.7, side=3, line=2, at=min(win_sites$pos[win_sites$win_id == win_id]) - .18 * (max(win_sites$pos[win_sites$win_id == win_id]) - min(win_sites$pos[win_sites$win_id == win_id])), adj=0)
  mtext(paste(round(sum(d_wins$d_plus[d_wins$win_id == win_id] > d_wins$d_plus) / length(d_wins$d_plus), 3), ' percentile', sep=''), 0, cex=.7, side=3, line=1, at=min(win_sites$pos[win_sites$win_id == win_id]) - .18 * (max(win_sites$pos[win_sites$win_id == win_id]) - min(win_sites$pos[win_sites$win_id == win_id])), adj=0)
  if (sum(genes$win_id == win_id) > 0)
    {
    rect(genes$start[genes$win_id == win_id], -1, genes$end[genes$win_id == win_id], 2, col=adjustcolor('gray', .2), border=NA)
    mtext(substr(genes$synonyms[genes$win_id == win_id], 1, 30), side = 3, line = 0:2, at = (sapply(genes$start[genes$win_id == win_id], max, min(win_sites$pos[win_sites$win_id == win_id])) + sapply(genes$end[genes$win_id == win_id], min, max(win_sites$pos[win_sites$win_id == win_id]))) / 2, cex=.5)
    }
  
  for (i in 1:nrow(sech_ssh_pops))
    {
    mtext(sech_ssh_pops$short_desc[i], side=2, at = i/nrow(sech_ssh_pops) - 1/nrow(sech_ssh_pops), line=5, col=sech_ssh_pops$col[i], las=2, cex=.5, adj=1)
    }

  #get sech alleles and adjusted allele distance for each sample
  par(mar=c(5.1, 10.1, 0, 1.1))
  alleles <- dbGetQuery(conn, paste("select a.pos, a.sample_id, a.num_der_alleles, abs(a2.num_der_alleles - a.num_der_alleles) num_der_alleles_adj
                                     from ", allele_table_sech, " a, ", allele_table_sech, " a2
                                     where a.win_id = '", win_id, "'
                                     and a.win_id = a2.win_id
                                     and a.pos = a2.pos
                                     and a.num_der_alleles != -2
                                     and a2.num_der_alleles != -2
                                     and a2.sample_id = 'SECH_3-sech_Anro_B3_TTAGGC_L001'
                                     union all
                                     select a.pos, a.sample_id, a.num_der_alleles, abs(a2.num_der_alleles - a.num_der_alleles) num_der_alleles_adj
                                     from ", allele_table_ssh, " a, ", allele_table_sech, " a2
                                     where a.win_id = '", win_id, "'
                                     and a.win_id = a2.win_id
                                     and a.pos = a2.pos
                                     and a.num_der_alleles != -2
                                     and a2.num_der_alleles != -2
                                     and a2.sample_id = 'SECH_3-sech_Anro_B3_TTAGGC_L001'
                                     order by a.pos, a.sample_id", sep=''))
  
  adj_allele_dist <- dbGetQuery(conn, paste("select a.sample_id, round(adj_allele_dist) + 1 adj_dist
                                             from ", allele_dist_table, " a, sample_pop_link l
                                             where win_id = '", win_id, "'
                                             and a.sample_id = l.sample_id
                                             and pop in (", pop_sql_str, ")
                                             order by l.pop, a.sample_id", sep=''))

  
  axis_cols <- merge(sech_ssh_pops, samples, by = 'pop')
  dist_cols <- matlab.like(101)
  
  #plot colored alleles for each sample. Alleles are polarized for a sech Anro sample (SECH_3-sech_Anro_B3_TTAGGC_L001) taking advantage of very low polymorphism.
  #the plots help identify which samples are polymorphic for the outlier windows and where haplotypic variation exists
  #blue: homozygous anro allele, red: homozygous alt allele, green heterozygous
  allele_cols <- adjustcolor(c('navy', 'forestgreen', 'maroon'), .6)
  geno_types = c('homo 1', 'het', 'homo 2')
  plot(1, type='n', xlim=range(win_sites$pos[win_sites$win_id == win_id]), ylim=c(0, nrow(samples)), xlab=paste('Pos on', d_wins$chrom[d_wins$win_id == win_id]), ylab='', main='', yaxt='n')
  for (i in 1:nrow(samples))
    {
    axis(2, at = i - 0.5, labels=samples$sample_id[i], las=2, cex.axis=.5, col.axis=axis_cols$col[i])
    points(alleles$pos[alleles$sample_id == samples$sample_id[i]], rep(i - 0.5, sum(alleles$sample_id == samples$sample_id[i])), pch=15, cex=.7, col=allele_cols[1 + alleles$num_der_alleles_adj[alleles$sample_id == samples$sample_id[i]]])
    rect(min(alleles$pos) - 3/25 * win_size , i - 1, min(alleles$pos) - 0.1/25 * win_size, i, col=dist_cols[adj_allele_dist$adj_dist[i]], border=NA)
    }
  
  for (i in 1:length(allele_cols))
    {text(min(alleles$pos) + (i-1) / 10 * (max(alleles$pos) - min(alleles$pos)), nrow(samples) + 1, geno_types[i], col=allele_cols[i])}
      
  abline(v=gwas_pos$pos[gwas_pos$win_id == win_id], col='black')

  #colored boxes are also plotted on the left side to indicate what proportion of the alleles differ from SECH_3-sech_Anro_B3_TTAGGC_L001 (blue low, red high)
  if (sum(genes$win_id == win_id) > 0)
    {
    rect(genes$start[genes$win_id == win_id], -nrow(samples), genes$end[genes$win_id == win_id], 1.5*nrow(samples), col=adjustcolor('gray', .2), border=NA)
    }  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  }  
layout(1)
dev.off()

