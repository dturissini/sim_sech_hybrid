library("RSQLite")  
library(colorRamps)

#process command line arguments
args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]


#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
pdf_file <- paste('sech_denis_dplus_win_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)

#load windows from db
d_stats_table <- paste('d_stat_win_', win_size, sep='')
wins <- dbGetQuery(conn, paste("select * 
                                from ", d_stats_table, "
                                where d_plus is not null", sep=''))

#get unique chromosomes
chroms <- sort(unique(wins$chrom))

#set threhsold for number of called sites per window
num_sites_cutoff <- as.numeric(win_size) / 50


pdf(pdf_file, height=8, width=10.5)
#hist of sites per window
hist(wins$num_sites, breaks=seq(0, max(wins$num_sites) + 100, 100), col='black', xlab='num sites', ylab='windows', main=c('Total sites per window', paste('window size =', win_size)))
abline(v=as.numeric(win_size) / 5, col='red')

#hist of D+
hist(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], breaks=100, col='black', xlab='D+', ylab='windows', main=c(paste('D+ for windows with > ', as.numeric(win_size) / 5, ' sites', sep=''), paste('window size =', win_size)))

#traces of D+ windows for each chrom
for (chrom in chroms)
  {
  win_filter <- wins$chrom == chrom & wins$num_sites > as.numeric(win_size) / 5
  plot((wins$start[win_filter] + wins$end[win_filter]) / 2, wins$d_plus[win_filter], type='l', ylim=c(min(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5]), max(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5])), xlab='pos', ylab='D+', main=paste('D+', 'window size =', win_size, chrom))
  abline(h=0, col='grey')
  }

#scatterplots of D+ vs site pattern differences
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], (wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5], pch=20, xlab='D+', ylab='ABBA - BABA', main=paste('D+ vs ABBA - BABA for window size =', win_size))
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], (wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5], pch=20, xlab='D+', ylab='BAAA - ABAA', main=paste('D+ vs BAAA - ABAA for window size =', win_size))


#abba - baba vs baaa - abaa
#colored by D+ value
d_plus_col_pal <- adjustcolor(matlab.like(101), .4)
d_plus_cols <- d_plus_col_pal[1 + round(100 * (wins$d_plus[wins$num_sites > as.numeric(win_size) / 5] - min(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5])) / (max(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5]) - min(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5])))]
plot((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5], (wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5], pch=20, col=d_plus_cols, xlab='ABBA - BABA', ylab='BAAA - ABAA', main=paste('ABBA - BABA vs BAAA - ABAA for window size =', win_size))

abline(v=0, col='gray')
abline(h=0, col='gray')


#plot custom legend
x_len <- max((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5]) - min((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5]) 
y_len <- max((wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5]) - min((wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5]) 

x_min <- min((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5])
x_max <- x_min + 0.05 * x_len
y_min <- min((wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5]) + 0.7 * y_len
y_max <- y_min + 0.25 * y_len
y_step <- (y_max - y_min) / length(d_plus_col_pal)

for (i in 1:length(d_plus_cols))
  {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=d_plus_col_pal[i], border=d_plus_col_pal[i])}
  
legend_d_plus <- round(seq(min(round(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], 2)), max(round(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], 2)), length.out=5),2)
text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_d_plus)), legend_d_plus, cex=.6)
text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'D+')


#colored by D+ > or < 0
d_plus_binary_cols <- rep(adjustcolor('black', .4), sum(wins$num_sites > as.numeric(win_size) / 5))
d_plus_binary_cols[wins$d_plus[wins$num_sites > as.numeric(win_size) / 5] < 0] <- adjustcolor('red', .4)
plot((wins$abba - wins$baba)[wins$num_sites > as.numeric(win_size) / 5], (wins$baaa - wins$abaa)[wins$num_sites > as.numeric(win_size) / 5], pch=20, col=d_plus_binary_cols, xlab='ABBA - BABA', ylab='BAAA - ABAA', main=paste('ABBA - BABA vs BAAA - ABAA for window size =', win_size))
abline(v=0, col='gray')
abline(h=0, col='gray')
legend('topleft', c('D+ >= 0', 'D+ < 0'), fill=c('black', 'red'), border=NA)




#D+ vs % abba and baa sites
per_abba <- 100 * (wins$abba / (wins$abba + wins$baba + wins$baaa + wins$abaa))[wins$num_sites > as.numeric(win_size) / 5]
per_baaa <- 100 * (wins$baaa / (wins$abba + wins$baba + wins$baaa + wins$abaa))[wins$num_sites > as.numeric(win_size) / 5]
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], per_abba[wins$num_sites > as.numeric(win_size) / 5], pch=20, xlab='D+', ylab='% ABBA sites', main=paste('D+ vs % ABBA sites for window size =', win_size))
plot(wins$d_plus[wins$num_sites > as.numeric(win_size) / 5], per_baaa[wins$num_sites > as.numeric(win_size) / 5], pch=20, xlab='D+', ylab='% BAAA sites', main=paste('D+ vs % BAAA sites for window size =', win_size))
dev.off()

