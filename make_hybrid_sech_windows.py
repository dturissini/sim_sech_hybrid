import argparse
import sqlite3
import os
import tempfile
from datetime import datetime

#python3 /proj/matutelb/projects/gwas/scripts/make_hybrid_sech_windows.py /proj/matutelb/projects/gwas/genotype_datasets/sech_oa/sim_sech_diff_sites.db


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("db_file")

    return parser.parse_args()


def main():
  args = parse_args()
  db_file = args.db_file  

  conn = sqlite3.connect(db_file)
    
  create_cur = conn.cursor()                            
  chrom_cur = conn.cursor()                            
  diff_cur = conn.cursor()                            
  hybrid_cur = conn.cursor()                            

  create_cur.execute("drop table if exists hybrid_sech_win")
  create_cur.execute("""create table hybrid_sech_win
                 (hsw_id int primary key,
                  win_size int,
                  chrom varchar(20),
                  start int,
                  end int,
                  hybrid_line varchar(60),
                  num_sech_sites int,
                  total_diff_sites int)""")
   
  create_cur.execute("create index idx_hsa_se on hybrid_sech_win(start, end)")
 
 

  win_sizes = [10000, 50000, 100000]
  
  hsw_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as o:
    for win_size in win_sizes:
      chrom_query = chrom_cur.execute("""select chrom, max(pos)
                                         from sim_sech_diff_sites
                                         group by chrom""")
      
      for (chrom, max_pos) in chrom_query:
        print(win_size, chrom, datetime.now())
        for win_start in range(1, max_pos, win_size):
          win_end = win_start + win_size - 1
          
          diff_site_query = diff_cur.execute(f"""select count(*)
                                                 from sim_sech_diff_sites
                                                 where chrom = '{chrom}'
                                                 and pos between {win_start} and {win_end}
                                                 and num_sech >= 20
                                                 and num_sim >= 60
                                                 and sim_ref_freq = 1
                                                 and sech_ref_freq = 0""")
          
          for (num_diff_sites, ) in diff_site_query:
            hybrid_sech_site_query = hybrid_cur.execute(f"""select hybrid_line, count(*)
                                                            from hybrid_sech_alleles h, sim_sech_diff_sites d
                                                            where h.chrom = '{chrom}'
                                                            and h.pos between {win_start} and {win_end}
                                                            and h.chrom = d.chrom
                                                            and h.pos = d.pos
                                                            and num_sech >= 20
                                                            and num_sim >= 60
                                                            and sim_ref_freq = 1
                                                            and sech_ref_freq = 0 
                                                            group by hybrid_line""")
            
            for (hybrid_line, num_sech_sites) in hybrid_sech_site_query:
              hsw_id += 1
              o.write(f"{hsw_id}\t{win_size}\t{chrom}\t{win_start}\t{win_end}\t{hybrid_line}\t{num_sech_sites}\t{num_diff_sites}\n")                  

    conn.close()
    o.flush()     
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {o.name} hybrid_sech_win" """)
  

if __name__ == '__main__':
  main()
