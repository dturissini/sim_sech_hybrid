import sys
import gzip
import sqlite3

pop_file = str(sys.argv[1])
db_file = str(sys.argv[2])

#create sample_pop table
conn = sqlite3.connect(db_file)  
cur = conn.cursor()    
cur.execute("drop table if exists lk_pop")
cur.execute("""create table lk_pop
               (ls_id int primary key,
                pop varchar(50),
                sample_id varchar(50))""")

cur.execute(f"create index idx_ls_pop on lk_pop(pop)")
cur.execute(f"create index idx_ls_sample_id on lk_pop(sample_id)")


with open(pop_file, 'r') as s:
  for line in s:
    line = line.strip()
    pop, pop, locations = line.split("\t")
    
    if locations == 'ALL':
      cur.execute(f"""insert into lk_pop
                      select null, '{pop}', sample_id
                      from sample_pop
                      where pop = '{pop}'""")      
    else:
      location_str = "'" + "', '".join(locations.split('~')) + "'"
      
      cur.execute(f"""insert into lk_pop
                      select null, '{pop}', sample_id
                      from sample_pop
                      where pop = '{pop}'
                      and location in ({location_str})""")

conn.commit()
conn.close()        