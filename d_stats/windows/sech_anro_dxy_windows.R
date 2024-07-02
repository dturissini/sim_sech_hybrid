library("RSQLite")  
library(colorRamps)


args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]



#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('sech_anro_dxy_win_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)

#load window and gwas data from db

dxy_wins_table <- paste('sech_anro_dxy_win_', win_size, sep='')
outlier_adj_allele_table <- paste('outlier_pi_sech_win_sech_adj_allele_dist_', win_size, sep='')
poly_win_table <- paste('poly_win_', win_size, sep='')


poly_wins <- dbGetQuery(conn, paste("select * 
                                     from ", poly_win_table, "
                                     where pi_sim is not null
                                     and pi_sech is not null
                                     and pi_ssh is not null", sep=''))



dxy_wins <- dbGetQuery(conn, paste("select s.sample_id, location, chrom, start, end, dxy_sech_anro, num_sites,
                                   case when location = 'Denis, Seychelles' then 'purple'
                                        when location = 'La Digue, Seychelles' then 'green'
                                        when location = 'Marianne, Seychelles' then 'darkblue'
                                        when location = 'Praslin, Seychelles' then 'lightblue'
                                   else 'gray' end col
                                   from sample_pop s, ", dxy_wins_table, " d, 
                                        (select sum(adj_allele_dist) cum_dist, sample_id
                                         from ", outlier_adj_allele_table, "
                                         group by sample_id) x
                                   where d.sample_id = s.sample_id
                                   and d.sample_id = x.sample_id
                                   and dxy_sech_anro is not null
                                   order by s.sample_id, chrom, start", sep=''))

win_avg_dxy <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist,
                                       case when location = 'Denis, Seychelles' then 'purple'
                                            when location = 'La Digue, Seychelles' then 'green'
                                            when location = 'Marianne, Seychelles' then 'darkblue'
                                            when location = 'Praslin, Seychelles' then 'lightblue'
                                       else 'gray' end col
                                       from sample_pop s, ", dxy_wins_table, " d, 
                                            (select sum(adj_allele_dist) cum_dist, sample_id
                                             from ", outlier_adj_allele_table, "
                                             group by sample_id) x
                                       where d.sample_id = s.sample_id
                                       and d.sample_id = x.sample_id
                                       group by s.sample_id, location", sep=''))


avg_dxy <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist,
                                   case when location = 'Denis, Seychelles' then 'purple'
                                        when location = 'La Digue, Seychelles' then 'green'
                                        when location = 'Marianne, Seychelles' then 'darkblue'
                                        when location = 'Praslin, Seychelles' then 'lightblue'
                                   else 'gray' end col
                                   from sample_pop s, ", dxy_wins_table, " d, 
                                        (select sum(adj_allele_dist) cum_dist, sample_id
                                         from ", outlier_adj_allele_table, "
                                         group by sample_id) x
                                   where d.sample_id = s.sample_id
                                   and d.sample_id = x.sample_id
                                   group by s.sample_id, location", sep=''))

avg_dxy_no_outlier_025_sech_pi_win <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist,
                                                              case when location = 'Denis, Seychelles' then 'purple'
                                                                   when location = 'La Digue, Seychelles' then 'green'
                                                                   when location = 'Marianne, Seychelles' then 'darkblue'
                                                                   when location = 'Praslin, Seychelles' then 'lightblue'
                                                              else 'gray' end col
                                                              from sample_pop s, ", dxy_wins_table, " d, 
                                                                   (select sum(adj_allele_dist) cum_dist, sample_id
                                                                    from ", outlier_adj_allele_table, "
                                                                    group by sample_id) x
                                                              where d.sample_id = s.sample_id
                                                              and d.sample_id = x.sample_id
                                                              and not exists (select 'x'
                                                                              from ", poly_win_table, " p
                                                                              where p.chrom = d.chrom
                                                                              and p.start = d.start
                                                                              and pi_sech > 0.025)
                                                              group by s.sample_id, location", sep=''))


avg_dxy_no_outlier_01_sech_pi_win <- dbGetQuery(conn, paste("select s.sample_id, location, sum(dxy_sech_anro) / count(*) avg_dxy, cum_dist,
                                                              case when location = 'Denis, Seychelles' then 'purple'
                                                                   when location = 'La Digue, Seychelles' then 'green'
                                                                   when location = 'Marianne, Seychelles' then 'darkblue'
                                                                   when location = 'Praslin, Seychelles' then 'lightblue'
                                                              else 'gray' end col
                                                              from sample_pop s, ", dxy_wins_table, " d, 
                                                                   (select sum(adj_allele_dist) cum_dist, sample_id
                                                                    from ", outlier_adj_allele_table, "
                                                                    group by sample_id) x
                                                              where d.sample_id = s.sample_id
                                                              and d.sample_id = x.sample_id
                                                              and not exists (select 'x'
                                                                              from ", poly_win_table, " p
                                                                              where p.chrom = d.chrom
                                                                              and p.start = d.start
                                                                              and pi_sech > 0.01)
                                                              group by s.sample_id, location", sep=''))


samples <- dbGetQuery(conn, paste("select sample_id, location
                                   from sample_pop s
                                   where sample_id in (select distinct sample_id from ", dxy_wins_table, ")
                                   order by location, sample_id", sep=''))


sample_loc <- data.frame(cbind(c('sech - Denis', 'sech - La Digue', 'sech - Marianne', 'sech - Praslin'),
                               c('Denis, Seychelles', 'La Digue, Seychelles', 'Marianne, Seychelles', 'Praslin, Seychelles'),
                               c('purple', 'green', 'darkblue', 'lightblue')))
names(sample_loc) <- c('display_loc', 'location', 'loc_col')


chroms <- sort(unique(dxy_wins$chrom))
num_sites_cutoff <- as.numeric(win_size) / 50


pdf(pdf_file, height=8, width=10.5)
plot(avg_dxy$avg_dxy, avg_dxy$cum_dist, pch=20, col=avg_dxy$col, xlab = 'Sech - Anro dxy', ylab = 'Anro allele distance', main = 'Sech dxy vs Anro allele distance')
legend('topleft', sample_loc$display_loc, fill=sample_loc$loc_col, border=NA)

plot(avg_dxy_no_outlier_025_sech_pi_win$avg_dxy, avg_dxy_no_outlier_025_sech_pi_win$cum_dist, pch=20, col=avg_dxy_no_outlier_025_sech_pi_win$col, xlab = 'Sech - Anro dxy', ylab = 'Anro allele distance', main = c('Sech dxy vs Anro allele distance for pi sech peaks', 'no outlier pi sech windows (pi sech > 0.025)'))
legend('topleft', sample_loc$display_loc, fill=sample_loc$loc_col, border=NA)

plot(avg_dxy_no_outlier_01_sech_pi_win$avg_dxy, avg_dxy_no_outlier_01_sech_pi_win$cum_dist, pch=20, col=avg_dxy_no_outlier_01_sech_pi_win$col, xlab = 'Sech - Anro dxy', ylab = 'Anro allele distance', main = c('Sech dxy vs Anro allele distance for pi sech peaks', 'no outlier pi sech windows (pi sech > 0.01)'))
legend('topleft', sample_loc$display_loc, fill=sample_loc$loc_col, border=NA)

for (chrom in chroms)
  {
  xlim_range <- c(0, max(poly_wins$end[poly_wins$chrom == chrom]))  
  
  par(mfrow=c(5,1), mar=c(0, 4.1, 4.1, 2.1))
  win_filter <- poly_wins$chrom == chrom & poly_wins$num_sites > num_sites_cutoff
  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$pi_sech[win_filter], xlim=xlim_range, type='l', col='red', xlab='', xaxt='n', ylab='Pi', main=chrom)


  y_max <- max(dxy_wins$dxy_sech_anro[dxy_wins$chrom == chrom  & dxy_wins$num_sites > num_sites_cutoff])  
  par(mar=c(0, 4.1, 0, 2.1))
  for (i in 1:nrow(sample_loc))
    {
    if (i == nrow(sample_loc))
      {par(mar=c(5.1, 4.1, 0, 2.1))}
    
    plot(1, type='n', xlim=xlim_range, ylim=c(0, y_max), xlab='', xaxt='n', ylab='Anro sech dxy', main='')
    text(1000000, .9* y_max, sample_loc$location[i])
    for (sample_id in samples$sample_id[samples$location == sample_loc$location[i]])
      {
      sample_win_filter <- dxy_wins$chrom == chrom & dxy_wins$num_sites > num_sites_cutoff & dxy_wins$sample_id == sample_id
    
      points((dxy_wins$start[sample_win_filter] + dxy_wins$end[sample_win_filter]) / 2, dxy_wins$dxy_sech_anro[sample_win_filter], type='l', col=dxy_wins$col[sample_win_filter])
      if (i == nrow(sample_loc))
        {axis(1, at=seq(0, max(dxy_wins$end[sample_win_filter]), 5000000))}
      }
    }
  
  par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1))
  }


par(mfrow=c(3,3))
for (sample_id in samples$sample_id)
  {
  hist(dxy_wins$dxy_sech_anro[dxy_wins$sample_id == sample_id], breaks=seq(0, max(dxy_wins$dxy_sech_anro[dxy_wins$sample_id == sample_id]) + .0005, .0005), col=dxy_wins$col[dxy_wins$sample_id == sample_id][1], border=NA, xlab='dxy - Anro', ylab='windows', main=c(sample_id, samples$location[samples$sample_id == sample_id]))
  }
par(mfrow=c(1,1))
dev.off()

