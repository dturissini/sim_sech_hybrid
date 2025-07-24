#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/assembly'
read_length_file <- paste(base_dir, '/sech_fastq_read_lengths.txt', sep='')

setwd(base_dir)

read_lens <- read.table(read_length_file, col.names=c('line', 'read_len', 'num_reads'))


pdf('sech_fastq_read_lengths.pdf', height=8, width=10.5)
for (line in sort(unique(read_lens$line)))
  {
  total_reads_line <- sum(read_lens$num_reads[read_lens$line == line])
  cum_seq_line_gb <- round(sum(read_lens$read_len[read_lens$line == line]) / 1e9, 2)
  
##  plot(1, type='n', xlim=c(0, max(read_lens$read_len)), ylim=c(0, max(read_lens$num_reads[read_lens$line == line])), xlab='Read length', ylab='Number of reads', main=c(line, paste('Total reads = ', total_reads_line, '; Cumulative seq length = ', cum_seq_line_gb, 'Gb', sep='')))
  plot(1, type='n', xlim=c(0, 10000), ylim=c(0, max(read_lens$num_reads[read_lens$line == line])), xlab='Read length', ylab='Number of reads', main=c(line, paste('Total reads = ', total_reads_line, '; Cumulative seq length = ', cum_seq_line_gb, 'Gb', sep='')))
  poly_xs <- c(read_lens$read_len[read_lens$line == line][1], read_lens$read_len[read_lens$line == line], read_lens$read_len[read_lens$line == line][sum(read_lens$line == line)], read_lens$read_len[read_lens$line == line][1])
  poly_ys <- c(0, read_lens$num_reads[read_lens$line == line], 0, 0)

  polygon(poly_xs, poly_ys, col='black') 
  }
dev.off()


cat('line', 'cum_len', 'num_reads', 'avg_len', 'per_gt_1kb', 'per_gt_5kb', 'per_gt_10kb',"\n")

for (line in sort(unique(read_lens$line)))
  {
  cum_len <- sum(read_lens$read_len[read_lens$line == line]) / 1e9
  num_reads <- sum(read_lens$num_reads[read_lens$line == line])
  avg_len <- round(sum(read_lens$read_len[read_lens$line == line] * read_lens$num_reads[read_lens$line == line]) / sum(read_lens$num_reads[read_lens$line == line]), 1)
  per_gt_1kb <- round(100 * sum(read_lens$num_reads[read_lens$line == line & read_lens$read_len > 1000]) / sum(read_lens$num_reads[read_lens$line == line]), 1)
  per_gt_5kb <- round(100 * sum(read_lens$num_reads[read_lens$line == line & read_lens$read_len > 5000]) / sum(read_lens$num_reads[read_lens$line == line]), 1)
  per_gt_10kb <- round(100 * sum(read_lens$num_reads[read_lens$line == line & read_lens$read_len > 10000]) / sum(read_lens$num_reads[read_lens$line == line]), 1)
  
  cat(line, cum_len, num_reads, avg_len, per_gt_1kb, per_gt_5kb, per_gt_10kb,"\n")
  }

line                          cum_len   num_reads avg_len per_gt_1kb per_gt_5kb per_gt_10kb 
DenisNF134                    0.2842109 9233658   700.4   14.6      1.6         0.3 
SECH_20-sech_PNF5_GTGGCC_L001 0.5476745 40789302  661.8   13.6      0.8         0.2 
SECH_21-sech_PNF3_GTTTCG_L001 0.6802232 77505052  602.3   12.3      0.6         0.1 
sech_Anro_114_female          1.874175  79258495  669.4   13        1           0.3 
sech_Denis_NF22               2.541594  80959491  728     13.6      1.3         0.5 
sech_Denis_NF45_female        2.656919  81310397  736.5   13.7      1.3         0.5 
sech_Denis_NF49               2.885583  82157947  762.7   13.9      1.5         0.6 
sech_Denis_NF69               3.178363  83803524  809.4   14.5      1.7         0.7 
sech_Denis_NF712              3.503808  84704243  852.6   14.9      1.9         0.9 
sech_Denis_papaya152          3.508355  84850156  854.2   14.9      1.9         0.9 
