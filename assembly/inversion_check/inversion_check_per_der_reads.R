library("RSQLite")  


#define file paths
base_dir <- '/work/users/d/t/dturissi/drosophila/ssh/assembly/inversion_check'
db_file <- paste(base_dir, '/', 'inversion_check_per_der.db', sep='')
pdf_file <- paste('inversion_check_per_der_reads.pdf', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)


sample_ids <- dbGetQuery(conn,"select sample_id
                               from inversion_check_per_der_reads
                               group by sample_id")



#plot derived allele freqs for each outlier window
pdf(pdf_file, height=8, width=10.5)
for (sample_id in sample_ids$sample_id)
  {
  anro_concord <- dbGetQuery(conn, paste("select inv_name, read_name, total_sites, 1.0 *num_anro_b3_alleles / total_sites anro_concordance
                                          from inversion_check_per_der_reads
                                          where sample_id = '", sample_id, "'
                                          and total_sites > 0", sep=''))

  hist(anro_concord$anro_concordance[substr(anro_concord$inv_name, 1, 4) != 'rand'], breaks = seq(0, 1.05, .05), col='black', xlab='Anro allele concordance', ylab='reads', main=c(sample_id, 'putative inversions'))
  plot(anro_concord$anro_concordance[substr(anro_concord$inv_name, 1, 4) != 'rand'], anro_concord$total_sites[substr(anro_concord$inv_name, 1, 4) != 'rand'], pch=20, col=adjustcolor('gray', .4), xlab='Anro allele concordance', ylab='total sites', main=c(sample_id, 'putative inversions'))
  
  if (sum(substr(anro_concord$inv_name, 1, 4) == 'rand') > 0)
    {
    hist(anro_concord$anro_concordance[substr(anro_concord$inv_name, 1, 4) == 'rand'], breaks = seq(0, 1.05, .05), col='black', xlab='Anro allele concordance', ylab='reads', main=c(sample_id, 'random regions'))
    plot(anro_concord$anro_concordance[substr(anro_concord$inv_name, 1, 4) == 'rand'], anro_concord$total_sites[substr(anro_concord$inv_name, 1, 4) == 'rand'], pch=20, col=adjustcolor('gray', .4), xlab='Anro allele concordance', ylab='total sites', main=c(sample_id, 'random regions'))
    }
    
  for (inv_name in sort(unique(anro_concord$inv_name)))
    {
    hist(anro_concord$anro_concordance[anro_concord$inv_name == inv_name], breaks = seq(0, 1.05, .05), col='black', xlab='Anro allele concordance', ylab='reads', main=c(sample_id, inv_name))   
    }
  }

dev.off()







