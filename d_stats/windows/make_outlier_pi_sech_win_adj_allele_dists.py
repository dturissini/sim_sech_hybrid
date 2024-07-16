# Import packages.
import sqlite3
import sys



def main():
  win_size = int(sys.argv[1])
  outlier_type = str(sys.argv[2])
  pop_str = str(sys.argv[3])
  db_file = str(sys.argv[4])
  
  #define database tables
  conn = sqlite3.connect(db_file)  
  win_allele_table_sech = "outlier_" + outlier_type + "_win_alleles_sech_" + str(win_size) + '_' + pop_str
  win_allele_table_ssh = "outlier_" + outlier_type + "_win_alleles_ssh_" + str(win_size) + '_' + pop_str
  adj_allele_dist_table = "outlier_" + outlier_type + "_win_sech_adj_allele_dist_" + str(win_size) + '_' + pop_str
  
  #create adj_allele_dist_table table
  cur = conn.cursor()    
  cur.execute(f"drop table if exists {adj_allele_dist_table}")
  cur.execute(f"""create table {adj_allele_dist_table}
                  (aad_id int primary key,
                   win_id varchar(30),
                   sample_id varchar(50),
                   adj_allele_dist float)""")

                   
  cur.execute(f"create index idx_add_win_id_{win_size}{outlier_type} on {adj_allele_dist_table}(win_id)")
  cur.execute(f"create index idx_add_sample_id_{win_size}{outlier_type} on {adj_allele_dist_table}(sample_id)")

  #check if ssh table exists
  ssh_table_exists = cur.execute(f"""select name 
                                     from sqlite_master 
                                     where type='table' 
                                     and name='{win_allele_table_ssh}'""").fetchone()
  
  #get adjusted allele distances for sech samples
  cur.execute(f"""insert into {adj_allele_dist_table}
                  select null, a.win_id, a.sample_id, 100.0 * (sum(abs(a2.num_der_alleles - a.num_der_alleles))) / count(*) / 2
                  from {win_allele_table_sech} a, {win_allele_table_sech} a2
                  where a.win_id = a2.win_id
                  and a.pos = a2.pos
                  and a.num_der_alleles != -2
                  and a2.num_der_alleles != -2
                  and a2.sample_id = 'SECH_3-sech_Anro_B3_TTAGGC_L001'
                  group by a.win_id, a.sample_id""")
  
  #get adjusted allele distances for ssh samples if an allele table exists for ssh
  if ssh_table_exists:
    cur.execute(f"""insert into {adj_allele_dist_table}
                    select null, a.win_id, a.sample_id, 100.0 * (sum(abs(a2.num_der_alleles - a.num_der_alleles))) / count(*) / 2
                    from {win_allele_table_ssh} a, {win_allele_table_sech} a2
                    where a.win_id = a2.win_id
                    and a.pos = a2.pos
                    and a.num_der_alleles != -2
                    and a2.num_der_alleles != -2
                    and a2.sample_id = 'SECH_3-sech_Anro_B3_TTAGGC_L001'
                    group by a.win_id, a.sample_id""")

  conn.commit()
   
if __name__ == '__main__':
  main()
