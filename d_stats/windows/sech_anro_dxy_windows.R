library("RSQLite")  

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]
pop_str <- args[2]



#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('sech_anro_dxy_win_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)


#define db tables
dxy_wins_table <- paste('sech_anro_dxy_win_', win_size, sep='')
outlier_adj_allele_table <- paste('outlier_pi_sech_win_sech_adj_allele_dist_', win_size, '_', pop_str, sep='')
poly_win_table <- paste('poly_win_', win_size, sep='')


#define pop str to use in db tables
pop_sql_str <- "'sechanro', 'sechden', 'sechlad', 'sechmari', 'sechpras', 'sechunk', 'sshlad', 'sshmahe', 'sshanro'"

#pops to exclude from pops due to low samples for sechunk and since sechanro is the basis of comparison for dxy
pop_exclude <- c('sechanro', 'sechunk')



#load window data from db
poly_wins <- dbGetQuery(conn, paste("select * 
                                     from ", poly_win_table, "
                                     where pop = 'sech'
                                     and win_id in (select distinct win_id 
                                                    from ", dxy_wins_table, ")", sep=''))



dxy_wins <- dbGetQuery(conn, paste("select s.sample_id, location, chrom, start, end, dxy_sech_anro, num_sites, col
                                   from sample_species s, ", dxy_wins_table, " d, pop_cols c, sample_pop_link l
                                   where d.sample_id = s.sample_id
                                   and dxy_sech_anro is not null
                                   and s.sample_id = l.sample_id
                                   and l.pop in (", pop_sql_str, ")
                                   and l.pop = c.pop
                                   order by s.sample_id, chrom, start", sep=''))




avg_dxy <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist, col
                                   from sample_species s, ", dxy_wins_table, " d, pop_cols c, sample_pop_link l,
                                        (select sum(adj_allele_dist) cum_dist, sample_id
                                         from ", outlier_adj_allele_table, "
                                         group by sample_id) x
                                   where d.sample_id = s.sample_id
                                   and d.sample_id = x.sample_id
                                   and s.sample_id = l.sample_id
                                   and l.pop in (", pop_sql_str, ")
                                   and l.pop = c.pop
                                   group by s.sample_id, location", sep=''))

avg_dxy_no_outlier_025_sech_pi_win <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist, col
                                                              from sample_species s, ", dxy_wins_table, " d, pop_cols c,  sample_pop_link l,
                                                                   (select sum(adj_allele_dist) cum_dist, sample_id
                                                                    from ", outlier_adj_allele_table, "
                                                                    group by sample_id) x
                                                              where d.sample_id = s.sample_id
                                                              and d.sample_id = x.sample_id
                                                              and not exists (select 'x'
                                                                              from ", poly_win_table, " p
                                                                              where p.chrom = d.chrom
                                                                              and p.start = d.start
                                                                              and pi > 0.025
                                                                              and pop = 'sech')
                                                              and s.sample_id = l.sample_id
                                                              and l.pop in (", pop_sql_str, ")
                                                              and l.pop = c.pop
                                                              group by s.sample_id, location", sep=''))


avg_dxy_no_outlier_01_sech_pi_win <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist, col
                                                              from sample_species s, ", dxy_wins_table, " d, pop_cols c,  sample_pop_link l,
                                                                   (select sum(adj_allele_dist) cum_dist, sample_id
                                                                    from ", outlier_adj_allele_table, "
                                                                    group by sample_id) x
                                                              where d.sample_id = s.sample_id
                                                              and d.sample_id = x.sample_id
                                                              and not exists (select 'x'
                                                                              from ", poly_win_table, " p
                                                                              where p.chrom = d.chrom
                                                                              and p.start = d.start
                                                                              and pi > 0.01
                                                                              and pop = 'sech')
                                                              and s.sample_id = l.sample_id
                                                              and l.pop in (", pop_sql_str, ")
                                                              and l.pop = c.pop
                                                              group by s.sample_id, location", sep=''))


#get samples, pops, and locations
samples <- dbGetQuery(conn, paste("select s.sample_id, location, l.pop
                                   from sample_species s, sample_pop_link l
                                   where s.sample_id = l.sample_id
                                   and l.pop in (", pop_sql_str, ")
                                   order by l.pop, s.sample_id", sep=''))


#get pops and plotting colors
sech_ssh_pops <- dbGetQuery(conn, paste("select pop, short_desc, col
                                         from pop_cols
                                         where pop in (", pop_sql_str, ")
                                         order by pop", sep=''))


#get unique chromosomes
chroms <- sort(unique(dxy_wins$chrom))

#set threshold of number of called sites per window
num_sites_cutoff <- as.numeric(win_size) / 50

#identify pi_sech 99% qunatile to identify outlier windows
pi_sech_threshold <- quantile(poly_wins$pi, .99, na.rm=T)



pdf(pdf_file, height=8, width=10.5)
#scatterplots comparing genome wide dxy with sech anro against the allele distance to a sech anro sample for outlier windows
#the idea is to see if samples that are more similar to sech from Anro for outlier windows are more similar to sech anro genome wide
plot(avg_dxy$avg_dxy, avg_dxy$cum_dist, pch=20, col=avg_dxy$col, xlab = 'Sech - Anro dxy', ylab = 'Anro allele distance', main = 'Sech dxy vs Anro allele distance')
legend('bottomright', sech_ssh_pops$short_desc, fill=sech_ssh_pops$col, border=NA)

plot(avg_dxy_no_outlier_025_sech_pi_win$avg_dxy, avg_dxy_no_outlier_025_sech_pi_win$cum_dist, pch=20, col=avg_dxy_no_outlier_025_sech_pi_win$col, xlab = 'Sech - Anro dxy', ylab = 'Anro allele distance', main = c('Sech dxy vs Anro allele distance for pi sech peaks', 'no outlier pi sech windows (pi sech > 0.025)'))
legend('bottomright', sech_ssh_pops$short_desc, fill=sech_ssh_pops$col, border=NA)

plot(avg_dxy_no_outlier_01_sech_pi_win$avg_dxy, avg_dxy_no_outlier_01_sech_pi_win$cum_dist, pch=20, col=avg_dxy_no_outlier_01_sech_pi_win$col, xlab = 'Sech - Anro dxy', ylab = 'Anro allele distance', main = c('Sech dxy vs Anro allele distance for pi sech peaks', 'no outlier pi sech windows (pi sech > 0.01)'))
legend('bottomright', sech_ssh_pops$short_desc, fill=sech_ssh_pops$col, border=NA)

#across each chromosome, plot pi sech and dxy for all pops compared with sech anro
for (chrom in chroms)
  {
  xlim_range <- c(0, max(poly_wins$end[poly_wins$chrom == chrom]))  
  
  par(mfrow=c(8,1), mar=c(0, 4.1, 4.1, 2.1))
  win_filter <- poly_wins$chrom == chrom & poly_wins$num_sites > num_sites_cutoff
  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$pi[win_filter], xlim=xlim_range, type='l', col='red', xlab='', xaxt='n', ylab='Pi', main=chrom)


  y_max <- max(dxy_wins$dxy_sech_anro[dxy_wins$chrom == chrom  & dxy_wins$num_sites > num_sites_cutoff])  
  par(mar=c(0, 4.1, 0, 2.1))
  for (i in which(!(sech_ssh_pops$pop %in% pop_exclude)))
    {
    if (i == max(which(!(sech_ssh_pops$pop %in% pop_exclude))))
      {par(mar=c(5.1, 4.1, 0, 2.1))}
    
    plot(1, type='n', xlim=xlim_range, ylim=c(0, y_max), xlab='', xaxt='n', ylab='Anro sech dxy', main='')
    text(1000000, .9* y_max, sech_ssh_pops$short_desc[i])
    for (sample_id in samples$sample_id[samples$pop == sech_ssh_pops$pop[i]])
      {
      sample_win_filter <- dxy_wins$chrom == chrom & dxy_wins$num_sites > num_sites_cutoff & dxy_wins$sample_id == sample_id
    
      points((dxy_wins$start[sample_win_filter] + dxy_wins$end[sample_win_filter]) / 2, dxy_wins$dxy_sech_anro[sample_win_filter], type='l', col=dxy_wins$col[sample_win_filter])
      if (i == nrow(sech_ssh_pops))
        {axis(1, at=seq(0, max(dxy_wins$end[sample_win_filter]), 5000000))}
      }
    
    rect(poly_wins$start[poly_wins$chrom == chrom & poly_wins$pi > pi_sech_threshold], -1, poly_wins$end[poly_wins$chrom == chrom & poly_wins$pi > pi_sech_threshold], 2, col=adjustcolor('gray', .2), border=NA)
    }
  
  par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1))
  }


#histograms of sech anro dxy windows for each sample
par(mfrow=c(3,3))
for (sample_id in samples$sample_id[samples$sample_id %in% sort(unique(dxy_wins$sample_id))])
  {
  hist(dxy_wins$dxy_sech_anro[dxy_wins$sample_id == sample_id], breaks=seq(0, max(dxy_wins$dxy_sech_anro[dxy_wins$sample_id == sample_id]) + .0005, .0005), col=dxy_wins$col[dxy_wins$sample_id == sample_id][1], border=NA, xlab='dxy - Anro', ylab='windows', main=c(sample_id, samples$location[samples$sample_id == sample_id]))
  }
par(mfrow=c(1,1))
dev.off()

