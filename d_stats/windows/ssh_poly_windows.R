library("RSQLite")  
library(fields)
library(colorRamps)

gridplot <- function(xvar, yvar, xbins, ybins, xbin_type = 'equal', ybin_type = 'equal', plot_type='normal', xlabel='', ylabel='', main_text='', z_lim=-1, bg='transparent', cell_counts=F, return_gp=F, print_legend=T, draw_axes=T, draw_plot=T)
  {  
  xbins_cut <- xbins
  ybins_cut <- ybins  
  xcut_right <- T
  ycut_right <- T
  
  if (xbin_type == 'equal')
    {            
    xbins_cut <- c(xbins[1] - diff(xbins)[1] / 2, xbins + c(diff(xbins)[1] / 2, diff(xbins) / 2))   
    }  else if (xbin_type == 'left')
    {
    xcut_right <- F  
    }  

  if (ybin_type == 'equal')
    {            
    ybins_cut <- c(ybins[1] - diff(ybins)[1] / 2, ybins + c(diff(ybins)[1] / 2, diff(ybins) / 2))   
    }  else if (ybin_type == 'left')
    {
    ycut_right <- F  
    }  
 
  tab <- table(cut(xvar, breaks=xbins_cut, right=xcut_right), cut(yvar, breaks=ybins_cut, right=ycut_right))             
  gp <- matrix(tab, ncol=length(xbins_cut) - 1, byrow=T)
                   
               
  if (plot_type == 'log')
    {
    plot_gp <- t(log(gp + 1, 10))  
#    legend_labels <- list(at = round(min(plot_gp),0):round(max(plot_gp),0), labels = 10**(round(min(plot_gp),0):round(max(plot_gp),0)))
    legend_labels <- list(at = round(min(plot_gp),0):round(max(plot_gp),0), labels = c(0, 10**(1:round(max(plot_gp),0))))
    
    if (print_legend == F)         
      {z_lim <- log(z_lim + 1, 10)}
    }  else if (plot_type == 'percent')
    {
    plot_gp <- t(gp / apply(gp,1,sum))
    legend_labels <- list(at = seq(0, 1, .1), labels = seq(0, 1, .1))
    }  else
    {
    plot_gp <- t(gp)    
    legend_labels <- NULL
    }
  
  if (draw_plot == T)
    {
    if (bg == 'black')  
      {par(bg='black', col.axis = 'white', col.lab='white', col.main = 'white', col.sub = 'white')}

    xaxes <- 's'
    yaxes <- 's'
    if (draw_axes == F)
      {
      xaxes <- 'n'
      yaxes <- 'n'
      }

#add row of negative counts to top of plot so only counts = 0 will appear as black.  Otherwise, low counts could be included in the black bin and appear to be 0.  
#likewise, ybins is modified to account for this extra row
#ylim is then set in the plots so the negative count row is not printed

#    plot_gp <- cbind(plot_gp, rep(-(max(plot_gp) - 32.5) / 64 + 1, nrow(plot_gp)))
    plot_gp <- cbind(plot_gp, rep(-(max(plot_gp) - 127.5) / 127, nrow(plot_gp)))
    ybins <- c(ybins, max(ybins) + diff(ybins)[length(ybins) - 1])

    
    if (print_legend == T)         
      {
      if (z_lim == -1) {
        image.plot(xbins, ybins, plot_gp, ylim=range(ybins_cut), axis.args = legend_labels, xaxt=xaxes, yaxt=yaxes, xlab = xlabel, ylab = ylabel, main=main_text, nlevel = 65, col=c('black', tim.colors()))
        } else {
        image.plot(xbins, ybins, plot_gp, ylim=range(ybins_cut), axis.args = legend_labels, xaxt=xaxes, yaxt=yaxes, xlab = xlabel, ylab = ylabel, main=main_text, nlevel = 65, zlim=z_lim, col=c('black', tim.colors()))          
        }
      } else if (z_lim == -1) {
      image(xbins, ybins, plot_gp, ylim=range(ybins_cut), xaxt=xaxes, yaxt=yaxes, xlab = xlabel, ylab = ylabel, main=main_text, nlevel = 65, col=c('black', tim.colors()))
      } else {
      image(xbins, ybins, plot_gp, ylim=range(ybins_cut), xaxt=xaxes, yaxt=yaxes, xlab = xlabel, ylab = ylabel, main=main_text, nlevel = 65, zlim=z_lim, col=c('black', tim.colors()))
      }

    
    if (cell_counts == T)
      {text(as.vector(sapply(xbins, rep, nrow(gp))), rep(ybins[1:length(ybins) - 1], ncol(gp)), gp)}            
     
    par(bg='transparent', col.axis = 'black', col.lab='black', col.main = 'black', col.sub = 'black')
    }

  if (return_gp == T)
    {plot_gp}     
  }


args <- commandArgs(trailingOnly = TRUE)
win_size <- args[1]

#define file paths
base_dir <- '/proj/matutelb/projects/drosophila/sim_sech_hybrid/introgression/d_stats/windows/'
db_file <- paste(base_dir, 'ssh_d_win.db', sep='')
gwas_snp_db_file <- '/proj/matutelb/projects/gwas/gwas_results/sech_oa_only_sim_pca/results/sech_oa_only_sim_pca_snp.db'
pdf_file <- paste('ssh_poly_win_', win_size, '.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, paste("attach database '", gwas_snp_db_file, "' as g", sep=''))

#load window and gwas data from db
d_stats_table <- paste('d_stat_win_', win_size, sep='')
poly_table <- paste('poly_win_', win_size, sep='')
d_wins <- dbGetQuery(conn, paste("select * 
                                  from ", d_stats_table, "
                                  where d_plus is not null", sep=''))

poly_wins <- dbGetQuery(conn, paste("select * 
                                     from ", poly_table, "
                                     where pi_sim is not null
                                     and pi_sech is not null
                                     and pi_ssh is not null", sep=''))

gwas_pos <- dbGetQuery(conn, paste("select g.chrom, pos, p, start, end, d_plus
                                    from gwas_snp g, ", d_stats_table, " w
                                    where g.chrom = w.chrom
                                    and pos between start and end
                                    and p < 1e-7", sep=''))

der_freq_ssh <- dbGetQuery(conn, paste("select pop, der_freq, num_sites
                                        from ssh_der_freq
                                        where der_freq > 0
                                        and der_freq < 1
                                        order by pop, der_freq", sep=''))

chroms <- sort(unique(d_wins$chrom))
num_sites_cutoff <- as.numeric(win_size) / 50


pdf(pdf_file, height=8, width=10.5)
#hist of sites per window
hist(d_wins$num_sites, breaks=seq(0, max(d_wins$num_sites) + 100, 100), col='black', xlab='num sites', ylab='windows', main=c('Total sites per window', paste('window size =', win_size)))
abline(v=num_sites_cutoff, col='red')

par(mfrow=c(3,1))
#pi
pi_breaks = seq(min(poly_wins$pi_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$pi_ssh[poly_wins$num_sites > num_sites_cutoff]), max(poly_wins$pi_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$pi_ssh[poly_wins$num_sites > num_sites_cutoff]) + .001, .001)
hist(poly_wins$pi_sim[poly_wins$num_sites > num_sites_cutoff], breaks=pi_breaks, col='black', xlab='pi', ylab='windows', main=c(paste('Pi simulans for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], breaks=pi_breaks, col='black', xlab='pi', ylab='windows', main=c(paste('Pi sechellia for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$pi_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=pi_breaks, col='black', xlab='pi', ylab='windows', main=c(paste('Pi sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))

#derived allele freq
der_freq_breaks = seq(min(poly_wins$der_freq_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_ssh[poly_wins$num_sites > num_sites_cutoff]), max(poly_wins$der_freq_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_ssh[poly_wins$num_sites > num_sites_cutoff]) + .001, .001)
hist(poly_wins$der_freq_sim[poly_wins$num_sites > num_sites_cutoff], breaks=der_freq_breaks, col='black', xlab='Derived allele freq', ylab='windows', main=c(paste('Derived freq simulans for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$der_freq_sech[poly_wins$num_sites > num_sites_cutoff], breaks=der_freq_breaks, col='black', xlab='Derived allele freq', ylab='windows', main=c(paste('Derived freq sechellia for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$der_freq_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=der_freq_breaks, col='black', xlab='Derived allele freq', ylab='windows', main=c(paste('Derived freq sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))


#derived allele freq diff
der_freq_diff_breaks = seq(min(poly_wins$der_freq_diff_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_diff_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_diff_sech_ssh[poly_wins$num_sites > num_sites_cutoff], na.rm=T), max(poly_wins$der_freq_diff_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_diff_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_diff_sech_ssh[poly_wins$num_sites > num_sites_cutoff], na.rm=T) + .001, .001)
hist(poly_wins$der_freq_diff_sim_sech[poly_wins$num_sites > num_sites_cutoff], breaks=der_freq_diff_breaks, col='black', xlab='Derived allele freq diff', ylab='windows', main=c(paste('Derived freq diff: simulans - sechellia for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$der_freq_diff_sim_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=der_freq_diff_breaks, col='black', xlab='Derived allele freq diff', ylab='windows', main=c(paste('Derived freq diff: simulans - sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$der_freq_diff_sech_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=der_freq_diff_breaks, col='black', xlab='Derived allele freq diff', ylab='windows', main=c(paste('Derived freq diff: sechellia - sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))


#dxy
dxy_breaks = seq(min(poly_wins$dxy_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sech_ssh[poly_wins$num_sites > num_sites_cutoff]), max(poly_wins$dxy_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sech_ssh[poly_wins$num_sites > num_sites_cutoff]) + .005, .005)
hist(poly_wins$dxy_sim_sech[poly_wins$num_sites > num_sites_cutoff], breaks=dxy_breaks, col='black', xlab='dxy', ylab='windows', main=c(paste('Dxy simulans - sechellia for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$dxy_sim_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=dxy_breaks, col='black', xlab='dxy', ylab='windows', main=c(paste('Dxy simulans - sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$dxy_sech_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=dxy_breaks, col='black', xlab='dxy', ylab='windows', main=c(paste('Dxy sechellia - sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))

#fst
fst_breaks = seq(min(poly_wins$fst_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$fst_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$fst_sech_ssh[poly_wins$num_sites > num_sites_cutoff]), max(poly_wins$fst_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$fst_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$fst_sech_ssh[poly_wins$num_sites > num_sites_cutoff]) + .005, .005)
hist(poly_wins$fst_sim_sech[poly_wins$num_sites > num_sites_cutoff], breaks=fst_breaks, col='black', xlab='Fst', ylab='windows', main=c(paste('Fst simulans - sechellia for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$fst_sim_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=fst_breaks, col='black', xlab='Fst', ylab='windows', main=c(paste('Fst simulans - sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
hist(poly_wins$fst_sech_ssh[poly_wins$num_sites > num_sites_cutoff], breaks=fst_breaks, col='black', xlab='Fst', ylab='windows', main=c(paste('Fst sechellia - sim_sech_hybrid for windows with > ', num_sites_cutoff, ' sites', sep=''), paste('window size =', win_size)))
par(mfrow=c(1,1))

#derived allele sfs
plot(der_freq_ssh$der_freq[der_freq_ssh$pop == 'sim'], der_freq_ssh$num_sites[der_freq_ssh$pop == 'sim'], type='l', col='blue', xlab='Derived freq', ylab='Sites', main='Per site derived allele freqs')
points(der_freq_ssh$der_freq[der_freq_ssh$pop == 'sech'], der_freq_ssh$num_sites[der_freq_ssh$pop == 'sech'], type='l', col='red')
points(der_freq_ssh$der_freq[der_freq_ssh$pop == 'ssh'], der_freq_ssh$num_sites[der_freq_ssh$pop == 'ssh'], type='l', col='gold')
legend('topright', c('sim', 'sech', 'ssh'), fill = c('blue', 'red', 'gold'), border=NA)


der_freq_range <- c(min(poly_wins$der_freq_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_ssh[poly_wins$num_sites > num_sites_cutoff]), max(poly_wins$der_freq_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_ssh[poly_wins$num_sites > num_sites_cutoff]))
plot(poly_wins$der_freq_sim[poly_wins$num_sites > num_sites_cutoff], poly_wins$der_freq_ssh[poly_wins$num_sites > num_sites_cutoff], pch=20, xlim=der_freq_range, ylim=der_freq_range, xlab='sim derived freq', ylab='sim_sech_hybrid derived freq', main='Derived allele freqs')
abline(0, 1, col='red')


bin_size <- .005
bins <- seq(floor(min(poly_wins$dxy_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sech_ssh[poly_wins$num_sites > num_sites_cutoff]) / bin_size) * bin_size, max(poly_wins$dxy_sim_sech[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sech_ssh[poly_wins$num_sites > num_sites_cutoff]) + bin_size, bin_size)
gridplot(poly_wins$dxy_sim_ssh[poly_wins$num_sites > num_sites_cutoff], poly_wins$dxy_sech_ssh[poly_wins$num_sites > num_sites_cutoff], bins, bins, xbin_type = 'left', ybin_type = 'left', plot_type='normal', xlabel='dxy: sim-ssh', ylabel='dxy: sech-ssh', main_text='dxy')
abline(0, 1, col='white')

#traces of D+ windows for each chrom
for (chrom in chroms)
  {
  xlim_range <- c(0, max(poly_wins$end[poly_wins$chrom == chrom]))  
  
  par(mfrow=c(7,1), mar=c(0, 4.1, 4.1, 2.1))
  gwas_p_range <- range(-log(gwas_pos$p, 10))
  plot(gwas_pos$pos[gwas_pos$chrom == chrom], -log(gwas_pos$p[gwas_pos$chrom == chrom], 10), pch=20, ylim=gwas_p_range, xlim=xlim_range, xlab='', xaxt='n', ylab='GWAS -log(Pi)', main=paste(chrom, "; win size = ", win_size, sep=''))
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)

  par(mar=c(0, 4.1, 0, 2.1))

  win_filter <- d_wins$chrom == chrom & d_wins$num_sites > num_sites_cutoff
  

  ab_range <- range(c(d_wins$abba[win_filter],  d_wins$baba[win_filter], d_wins$baaa[win_filter], d_wins$abaa[win_filter]))
  plot((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$abba[win_filter], col='darkblue', type='l', ylim=ab_range, xlim=xlim_range, xlab='', xaxt='n', ylab='AB site patterns', main='')
  points((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$baba[win_filter], col='lightblue', type='l')
  points((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$baaa[win_filter], col='purple', type='l')
  points((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$abaa[win_filter], col='pink', type='l')
  legend('topleft', c('ABBA', 'BABA', 'BAAA', 'ABAA'), fill = c('darkblue', 'lightblue', 'purple', 'pink'), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-9999, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(9999, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)



  plot((d_wins$start[win_filter] + d_wins$end[win_filter]) / 2, d_wins$d_plus[win_filter], type='l', ylim=c(min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]), max(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff])), xlim=xlim_range, xlab='', xaxt='n', ylab='D+', main='')
  abline(h=0, col='grey')
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)


  win_filter <- poly_wins$chrom == chrom & poly_wins$num_sites > num_sites_cutoff

  pi_range <- range(poly_wins$pi_sim[win_filter], poly_wins$pi_sech[win_filter], poly_wins$pi_ssh[win_filter])
  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$pi_sim[win_filter], xlim=xlim_range, type='l', col='blue', ylim=pi_range, xlab='', xaxt='n', ylab='Pi', main='')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$pi_sech[win_filter], type='l', col='red')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$pi_ssh[win_filter], type='l', col='gold')
  legend('topleft', c('sim', 'sech', 'ssh'), fill = c('blue', 'red', 'gold'), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)

  dxy_range <- range(poly_wins$dxy_sim_sech[win_filter], poly_wins$dxy_sim_ssh[win_filter], poly_wins$dxy_sech_ssh[win_filter])
  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$dxy_sim_sech[win_filter], xlim=xlim_range, type='l', col='purple', ylim=dxy_range, xlab='', xaxt='n', ylab='Dxy', main='')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$dxy_sim_ssh[win_filter], type='l', col='green')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$dxy_sech_ssh[win_filter], type='l', col='orange')
  legend('topleft', c('sim-sech', 'sim-ssh', 'sech-ssh'), fill = c('purple', 'green', 'orange'), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  
  fst_range <- range(poly_wins$fst_sim_sech[win_filter], poly_wins$fst_sim_ssh[win_filter], poly_wins$fst_sech_ssh[win_filter])
  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$fst_sim_sech[win_filter], xlim=xlim_range, type='l', col='purple', ylim=fst_range, xlab='Pos', xaxt='n', ylab='Fst', main='')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$fst_sim_ssh[win_filter], type='l', col='green')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$fst_sech_ssh[win_filter], type='l', col='orange')
  legend('topleft', c('sim-sech', 'sim-ssh', 'sech-ssh'), fill = c('purple', 'green', 'orange'), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)

  par(mar=c(5.1, 4.1, 0, 2.1))
  der_freq_range <- range(poly_wins$der_freq_sim[win_filter], poly_wins$der_freq_sech[win_filter], poly_wins$der_freq_ssh[win_filter])
  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$der_freq_sim[win_filter], xlim=xlim_range, type='l', col='blue', ylim=der_freq_range, xlab='', xaxt='n', ylab='Der freq', main='')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$der_freq_sech[win_filter], type='l', col='red')
  points((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, poly_wins$der_freq_ssh[win_filter], type='l', col='gold')
  legend('topleft', c('sim', 'sech', 'ssh'), fill = c('blue', 'red', 'gold'), border=NA)
  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  
  
##  T_sim_sech <- -log(1 - poly_wins$fst_sim_sech[win_filter], 10)
##  T_sim_ssh <- -log(1 - poly_wins$fst_sim_ssh[win_filter], 10)
##  T_sech_ssh <- -log(1 - poly_wins$fst_sech_ssh[win_filter], 10)
##  
##  PBS <- (T_sim_ssh + T_sech_ssh - T_sim_sech) / 2
##  
##  PBS_range <- range(PBS)
##  plot((poly_wins$start[win_filter] + poly_wins$end[win_filter]) / 2, PBS, xlim=xlim_range, type='l', ylim=PBS_range, xlab='Pos', ylab='PBS(ssh)', main='')
##  rect(gwas_pos$start[gwas_pos$chrom == chrom], rep(-99, sum(gwas_pos$chrom == chrom)), gwas_pos$end[gwas_pos$chrom == chrom], rep(99, sum(gwas_pos$chrom == chrom)), col=adjustcolor('gray', .2), border=NA)
  }
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))


plot(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], d_wins$d_plus[d_wins$num_sites > num_sites_cutoff], pch=20, col=adjustcolor('gray', .6), xlab='pi sech', ylab='D+', main='pi sech vs D+')

d_plus_col_pal <- adjustcolor(matlab.like(101), .4)
d_plus_cols <- d_plus_col_pal[1 + round(100 * (d_wins$d_plus[d_wins$num_sites > num_sites_cutoff] - min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff])) / (max(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) - min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff])))]
plot(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], d_wins$baaa[d_wins$num_sites > num_sites_cutoff] - d_wins$abaa[d_wins$num_sites > num_sites_cutoff], pch=20, col=d_plus_cols, xlab='pi sech', ylab='BAAA-ABAA', main='pi sech vs BAAA - ABAA')

#plot custom legend
x_len <- max(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff]) - min(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff]) 
y_len <- max((d_wins$baaa - d_wins$abaa)[d_wins$num_sites > num_sites_cutoff]) - min((d_wins$baaa - d_wins$abaa)[d_wins$num_sites > num_sites_cutoff]) 

x_min <- min(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff]) + .9 * x_len
x_max <- x_min + 0.05 * x_len
y_min <- min((d_wins$baaa - d_wins$abaa)[d_wins$num_sites > num_sites_cutoff]) + 0.1 * y_len
y_max <- y_min + 0.25 * y_len
y_step <- (y_max - y_min) / length(d_plus_col_pal)

for (i in 1:length(d_plus_cols))
  {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=d_plus_col_pal[i], border=d_plus_col_pal[i])}
  
legend_d_plus <- round(seq(min(round(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff], 2)), max(round(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff], 2)), length.out=5),2)
text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_d_plus)), legend_d_plus, cex=.6)
text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'D+')



sech_pi_col_pal <- adjustcolor(matlab.like(101), .4)
sech_pi_cols <- sech_pi_col_pal[1 + round(100 * (poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff] - min(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff])) / (max(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff]) - min(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff])))]
cexes <- 3 * poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff] / max(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff])
plot(d_wins$baaa[d_wins$num_sites > num_sites_cutoff] - d_wins$abaa[d_wins$num_sites > num_sites_cutoff], d_wins$d_plus[d_wins$num_sites > num_sites_cutoff], pch=20, cex=cexes, col=sech_pi_cols, xlab='BAAA-ABAA', ylab='D+', main='BAAA - ABAA vs D+')

#plot custom legend
x_len <- max((d_wins$baaa - d_wins$abaa)[d_wins$num_sites > num_sites_cutoff]) - min((d_wins$baaa - d_wins$abaa)[d_wins$num_sites > num_sites_cutoff]) 
y_len <- max(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) - min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) 

x_min <- min((d_wins$baaa - d_wins$abaa)[d_wins$num_sites > num_sites_cutoff]) + 0.9 * x_len
x_max <- x_min + 0.05 * x_len
y_min <- min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) + .05 * y_len
y_max <- y_min + 0.25 * y_len
y_step <- (y_max - y_min) / length(sech_pi_col_pal)

for (i in 1:length(sech_pi_cols))
  {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=sech_pi_col_pal[i], border=sech_pi_col_pal[i])}
  
legend_pi_sech <- round(seq(min(round(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], 2)), max(round(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], 2)), length.out=5),2)
text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_pi_sech)), legend_pi_sech, cex=.6)
text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'pi sech')


cexes <- 3 * poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff] / max(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff])
plot(d_wins$abba[d_wins$num_sites > num_sites_cutoff] - d_wins$baba[d_wins$num_sites > num_sites_cutoff], d_wins$d_plus[d_wins$num_sites > num_sites_cutoff], pch=20, cex=cexes, col=sech_pi_cols, xlab='ABBA-BABA', ylab='D+', main='ABBA - BABA vs D+')

#plot custom legend
x_len <- max((d_wins$abba - d_wins$baba)[d_wins$num_sites > num_sites_cutoff]) - min((d_wins$abba - d_wins$baba)[d_wins$num_sites > num_sites_cutoff]) 
y_len <- max(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) - min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) 

x_min <- min((d_wins$abba - d_wins$baba)[d_wins$num_sites > num_sites_cutoff]) + 0.9 * x_len
x_max <- x_min + 0.05 * x_len
y_min <- min(d_wins$d_plus[d_wins$num_sites > num_sites_cutoff]) + .05 * y_len
y_max <- y_min + 0.25 * y_len
y_step <- (y_max - y_min) / length(sech_pi_col_pal)

for (i in 1:length(sech_pi_cols))
  {rect(x_min, y_min + y_step * (i - 1), x_max, y_min + y_step * i, col=sech_pi_col_pal[i], border=sech_pi_col_pal[i])}
  
legend_pi_sech <- round(seq(min(round(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], 2)), max(round(poly_wins$pi_sech[poly_wins$num_sites > num_sites_cutoff], 2)), length.out=5),2)
text(x_max + 0.02 * x_len, seq(y_min, y_max, length.out=length(legend_pi_sech)), legend_pi_sech, cex=.6)
text(mean(c(x_min, x_max)), y_max + 0.03 * y_len, 'pi sech')
dev.off()

