library("RSQLite")  
library(colorRamps)
setwd('/work/users/d/t/dturissi/drosophila/ssh/introgression/pca')


pca_plots <- function (pca_name, pops_df)
  {
  pca_file <- paste(pca_name, '.eigenvec', sep='')
  pdf_file <- paste(pca_name, '.pdf', sep='')
  
  pca_raw <- read.table(pca_file, header=T, comment.char='')  
  pca <- merge(pca_raw, pops_df, by.x='IID', by.y= 'sample_id', all.x=T)
    
  pops <- unique(pops_df[, 2:3])  

  
  pdf(pdf_file, height=8, width=10.5)
  for (i in 1:9)
    {
    plot(pca[, i + 2], pca[, i + 3],  pch=20, col=pca$col, xlab=paste('PC', i, sep=''), ylab=paste('PC', i + 1, sep=''),  main='')
    legend('bottomright', pops$pop, fill=pops$col, border=NA)
    }
  dev.off()
  }


db_file <- '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
conn <- dbConnect(dbDriver("SQLite"), db_file)

#all
all_pops <- dbGetQuery(conn, "select sample_id, l.pop, col
                              from lk_pop l, pop_cols c
                              where l.pop in ('sim', 'sech', 'ssh', 'mel')
                              and l.pop = c.pop")
                              
pca_plots('sim_sech_outgroup_mel_pca', all_pops)



#sech pops
sech_pops <- dbGetQuery(conn, "select sample_id, l.pop, col
                              from lk_pop l, pop_cols c
                              where l.pop in ('sechden', 'sechpras', 'sechanro', 'sechlad', 'sechunk', 'sechmari')
                              and l.pop = c.pop")

pca_plots('sech_pca', sech_pops)


#sechbase pops
sechbase_pops <- dbGetQuery(conn, "select sample_id, l.pop, col
                                   from lk_pop l, pop_cols c
                                   where l.pop in ('sechanro', 'sechlad', 'sechmari')
                                   and l.pop = c.pop")

pca_plots('sechbase_pca', sechbase_pops)


#sech_ssh
sech_ssh_pops <- dbGetQuery(conn, "select sample_id, l.pop, col
                                   from lk_pop l, pop_cols c
                                   where l.pop in ('sshlad', 'sshmahe', 'sshanro', 'sechden', 'sechpras', 'sechanro', 'sechlad', 'sechunk', 'sechmari')
                                   and l.pop = c.pop")

pca_plots('sech_ssh_pca', sech_ssh_pops)
