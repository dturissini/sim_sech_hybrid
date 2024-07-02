# Import packages.
import sqlite3
import sys



def main():
  win_size = int(sys.argv[1])
  outlier_type = str(sys.argv[2])
  db_file = str(sys.argv[3])
  
  conn = sqlite3.connect(db_file)  
  allele_table = "outlier_" + outlier_type + "_win_sechellia_alleles_" + str(win_size)
  adj_allele_dist_table = "outlier_" + outlier_type + "_win_sech_adj_allele_dist_" + str(win_size)

  cur = conn.cursor()    
  cur.execute(f"drop table if exists {adj_allele_dist_table}")
  cur.execute(f"""create table {adj_allele_dist_table}
                  (aad_id int primary key,
                   dsw_id int,
                   sample_id varchar(50),
                   adj_allele_dist float)""")

                   
  cur.execute(f"create index idx_add_dsw_{win_size}{outlier_type} on {adj_allele_dist_table}(dsw_id)")
  cur.execute(f"create index idx_add_sample_id_{win_size}{outlier_type} on {adj_allele_dist_table}(sample_id)")

  cur.execute(f"""insert into {adj_allele_dist_table}
                  select null, a.dsw_id, a.sample_id, 100.0 * (sum(abs(a2.num_der_alleles - a.num_der_alleles))) / count(*) / 2
                  from {allele_table} a, {allele_table} a2
                  where a.dsw_id = a2.dsw_id
                  and a.pos = a2.pos
                  and a.num_der_alleles != -2
                  and a2.num_der_alleles != -2
                  and a2.sample_id = 'SECH_3-sech_Anro_B3_TTAGGC_L001'
                  group by a.dsw_id, a.sample_id""")

  conn.commit()
   
if __name__ == '__main__':
  main()
