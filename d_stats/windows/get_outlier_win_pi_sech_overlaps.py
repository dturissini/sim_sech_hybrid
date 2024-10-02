# Import packages.
import numpy as np
import numcodecs
import pandas as pd
import sys
import sqlite3


#create outlier overlap table
def create_windows_table(conn, win_size, a, b, overlap_table):
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {overlap_table}")
  cur.execute(f"""create table {overlap_table}
                  (odwps_id int primary key,
                   win_id varchar(30),
                   quantile decimal(4,3),
                   d_plus_outlier_a int,
                   d_plus_outlier_b int)""")

                   
  cur.execute(f"create index idx_odwps_win_id_{win_size}{a}{b} on {overlap_table}(win_id)")

  cur.close()
  


def main():
  win_size = int(sys.argv[1])
  pop_str_a = str(sys.argv[2])
  pop_str_b = str(sys.argv[3])
  db_file = str(sys.argv[4])
  
  
  
  #establish db connection and create d_stats_win table
  conn = sqlite3.connect(db_file)  

  pop_str_a = 'sim_ssh_sech_mel'
  win_size =50000
  pop_str_b = 'sechpras_sechbase_sim_mel'
  db_file = '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
  
  d_win_table_a = "d_stat_win_" + str(win_size) + '_' + pop_str_a
  d_win_table_b = "d_stat_win_" + str(win_size) + '_' + pop_str_b
  poly_win_table = "poly_win_" + str(win_size)
  overlap_table = "outlier_pi_sech_overlap_" + str(win_size) + '_' + pop_str_a + '_' + pop_str_b
  
  create_windows_table(conn, win_size, pop_str_a, pop_str_b, overlap_table)  
    
  quantiles = [.99, .98, .975, .95, .9]
  
  
  d_plus_a = pd.read_sql(f"""select win_id, d_plus
                             from {d_win_table_a}
                             where num_sites > 1000""", conn)
  
  d_plus_b = pd.read_sql(f"""select win_id, d_plus
                             from {d_win_table_b}
                             where num_sites > 1000""", conn)
  
  pi_sech = pd.read_sql(f"""select win_id, pi
                            from {poly_win_table}
                            where num_sites > 1000
                            and pop = 'sech'""", conn)
  

  for quantile in quantiles:
    pi_quantile = np.quantile(pi_sech['pi'], quantile)
    d_plus_quantile_a = np.quantile(d_plus_a['d_plus'], quantile)
    d_plus_quantile_b = np.quantile(d_plus_b['d_plus'], quantile)
    
    for win_id in pi_sech['win_id'][pi_sech['pi'] > pi_quantile]:
      overlap_a = sum([win_id in d_plus_a['win_id'][d_plus_a['d_plus'] > d_plus_quantile_a].values])
      overlap_b = sum([win_id in d_plus_a['win_id'][d_plus_b['d_plus'] > d_plus_quantile_b].values])
      
      conn.execute(f"""insert into {overlap_table}
                       values
                       (null, '{win_id}', {quantile}, {overlap_a}, {overlap_b})""")
      
  conn.commit()
  conn.close()
  

if __name__ == '__main__':
  main()
