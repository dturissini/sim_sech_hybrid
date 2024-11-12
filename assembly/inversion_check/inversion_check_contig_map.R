library("RSQLite")  

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- as.character(args[1])
sample_id <- args[2]
seq_type <- args[3]

#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/assembly/inversion_check'
pdf_file <- paste('inversion_check_', sample_id, '_', seq_type, '.pdf', sep='')

db_file <- 'minimap_results.db'
ssh_db_file <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'


setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
ssh_conn <- dbConnect(dbDriver("SQLite"), ssh_db_file)


possible_inversions <- read.table('denis_possible_inversions.txt', header=T, sep="\t")

#define db table
win_site_table <- paste('outlier_pi_sech_win_neighbor_sites_', win_size, '_sim_ssh_sech_mel', sep='')

contig_overlaps <- dbGetQuery(conn, paste("select inv_name, contig_name, contig_len, contig_start, contig_end, chrom, start, end
                                           from minimap_results
                                           where sample_id = '", sample_id, "'
                                           and seq_type = '", seq_type, "'", sep=''))


pdf(pdf_file, height=8, width=10.5)
for (i in which(substr(possible_inversions$name, 1, 4) != 'rand'))
  {
  layout(matrix(c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2), nrow=4, byrow=T))
  par(mar=c(0, 4.1, 4.1, 1.1))
  
  win_sites <- dbGetQuery(ssh_conn, paste("select *
                                           from ", win_site_table, "
                                           where der_alleles > 0
                                           and pop = 'sech'
                                           and chrom = '", possible_inversions$chrom[i], "'
                                           and pos between ", possible_inversions$start[i], " and ", possible_inversions$end[i], "
                                           order by pos", sep=''))


  plot(win_sites$pos, win_sites$der_alleles / win_sites$total_alleles, pch=20, col='red', xaxt='n', xlab='', ylim=c(0,1), ylab='Derived freq', main=possible_inversions$name[i])
  par(mar=c(5.1, 4.1, 0, 1.1))
  
  inv_overlaps <- subset(contig_overlaps, inv_name == possible_inversions$name[i])
  inv_overlaps_contigs <- aggregate(start ~ contig_name, inv_overlaps, min)
  inv_overlaps_contigs <- inv_overlaps_contigs[order(inv_overlaps_contigs$start), ]
  
  
  plot(1, type='n', xlim=c(min(win_sites$pos), max(win_sites$pos)), ylim=c(1, nrow(inv_overlaps)), xlab=possible_inversions$chrom[i], ylab='', yaxt='n', main='')
  for (j in nrow(inv_overlaps_contigs):1)
    {
    line_col <- 'black'    
    if (sum(inv_overlaps$contig_name == inv_overlaps_contigs$contig_name[j]) > 1)
      {line_col <- 'red'}
    
    for (row_id in which(inv_overlaps$contig_name == inv_overlaps_contigs$contig_name[j]))
      {
      segments(inv_overlaps$start[row_id], j, inv_overlaps$end[row_id], j, col=line_col)
      mtext(inv_overlaps$contig_name[row_id], side=2, at=j, cex=.4, las=1)
      }
    }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  layout(1)
  }  
dev.off()

