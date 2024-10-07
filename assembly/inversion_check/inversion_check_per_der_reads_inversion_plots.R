library("RSQLite")  
library(colorRamps)


base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/assembly/inversion_check'
db_file <- paste(base_dir, '/', 'inversion_check_per_der.db', sep='')
d_win_db_file <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
win_size <- 50000
possible_inversions <- read.table('denis_possible_inversions.txt', header=T, sep="\t")


win_site_table <- paste('outlier_pi_sech_win_neighbor_sites_', win_size, '_sim_ssh_sech_mel', sep='')


setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
d_win_conn <- dbConnect(dbDriver("SQLite"), d_win_db_file)


sample_ids <- dbGetQuery(conn,"select sample_id
                               from inversion_check_per_der_reads
                               group by sample_id")


concord_pal <- adjustcolor(matlab.like(101), .4)

for (sample_id in sample_ids$sample_id)
  {
  pdf_file <- paste('inversion_check_per_der_reads_inversion_plots_', sample_id, '.pdf', sep='')
  read_file <- paste(base_dir, '/inversion_overlaps_trimmedReads/', sample_id, '_inv_overlap_trimmedReads.txt', sep='')
  
  read_overlaps <- read.table(read_file, header=F, sep="\t")
  names(read_overlaps) <- c('inv_name', 'read_name', 'read_len', 'read_start', 'read_end', 'chrom', 'start', 'end')
  
  pdf(pdf_file, height=8, width=10.5)
  for (i in 1:nrow(possible_inversions))
    {
    layout(matrix(c(1, 2, 3, 4, 4, 4), nrow=6, byrow=T))
    par(mar=c(0, 4.1, 4.1, 1.1))
    
    win_sites <- dbGetQuery(d_win_conn, paste("select *
                                         from ", win_site_table, "
                                         where der_alleles > 0
                                         and pop = 'sech'
                                         and chrom = '", possible_inversions$chrom[i], "'
                                         and pos between ", possible_inversions$start[i], " and ", possible_inversions$end[i], "
                                         order by pos", sep=''))
  
  
    plot(win_sites$pos, win_sites$der_alleles / win_sites$total_alleles, pch=20, col='red', xaxt='n', xlab='', xlim=c(possible_inversions$start[i], possible_inversions$end[i]), ylim=c(0,1), ylab='Derived freq', main=possible_inversions$name[i])

    par(mar=c(0, 4.1, 0, 1.1))
    read_alleles <- dbGetQuery(conn, paste("select pos, 1.0 * sum(der_site) / count(*) der_freq, 1.0 * sum(anro_site) / count(*) anro_freq
                                                  from inversion_check_read_alleles
                                                  where sample_id = '", sample_id, "'
                                                  and inv_name = '", possible_inversions$name[i], "'
                                                  group by pos", sep=''))

    win_size <- 2000
    starts <- seq(possible_inversions$start[i], possible_inversions$end[i] + win_size, win_size)
    avg_win_der_freqs <- c()
    avg_win_anro_freqs <- c()
    for (start in starts)
      {
      win_filter <- read_alleles$pos >= start & read_alleles$pos < start + win_size
      avg_win_der_freqs <- c(avg_win_der_freqs, sum(read_alleles$der_freq[win_filter]) / sum(win_filter))
      avg_win_anro_freqs <- c(avg_win_anro_freqs, sum(read_alleles$anro_freq[win_filter]) / sum(win_filter))
      }
 
    plot(read_alleles$pos, read_alleles$der_freq, pch=20, xaxt='n', xlab='', xlim=c(possible_inversions$start[i], possible_inversions$end[i]), ylim=c(0,1), ylab='Nanopore derived freq')
    points(starts + win_size / 2, avg_win_der_freqs, type='l', col='blue')
    
    
    plot(read_alleles$pos, read_alleles$anro_freq, pch=20, xaxt='n', xlab='', xlim=c(possible_inversions$start[i], possible_inversions$end[i]), ylim=c(0,1), ylab='Nanopore Anro freq')
    points(starts + win_size / 2, avg_win_anro_freqs, type='l', col='red')
    
    par(mar=c(5.1, 4.1, 0, 1.1))
    
    inv_overlaps <- subset(read_overlaps, inv_name == possible_inversions$name[i])
    inv_overlaps_reads <- aggregate(start ~ read_name, inv_overlaps, min)
    inv_overlaps_reads <- inv_overlaps_reads[order(inv_overlaps_reads$start), ]
    
    
    plot(1, type='n', ylim=c(1, nrow(inv_overlaps)), xlab=possible_inversions$chrom[i], xlim=c(possible_inversions$start[i], possible_inversions$end[i]), ylab='', yaxt='n', main='')
    for (j in nrow(inv_overlaps_reads):1)
      {
      anro_concord <- dbGetQuery(conn, paste("select inv_name, read_name, total_sites, 1.0 * num_anro_b3_alleles / total_sites anro_concordance
                                              from inversion_check_per_der_reads
                                              where sample_id = '", sample_id, "'
                                              and inv_name = '", possible_inversions$name[i], "'
                                              and read_name = '", inv_overlaps_reads$read_name[j], "'
                                              and total_sites > 0", sep=''))
      
      col_i <- round(100 * anro_concord$anro_concordance)
      
#      cat(j, inv_overlaps_reads$read_name[j], col_i, anro_concord$anro_concordance, "\n")
      
      for (row_id in which(inv_overlaps$read_name == inv_overlaps_reads$read_name[j]))
        {
        segments(inv_overlaps$start[row_id], j, inv_overlaps$end[row_id], j, col=concord_pal[col_i])
        mtext(inv_overlaps$read_name[row_id], side=2, at=j, cex=.4, las=1)
        }
      }
    
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    layout(1)
    }  
  
  dev.off()
  }


