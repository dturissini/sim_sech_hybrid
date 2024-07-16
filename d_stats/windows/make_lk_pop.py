import sys
import gzip
import sqlite3

pop_file = str(sys.argv[1])
db_file = str(sys.argv[2])

#make sample_pop_link table
conn = sqlite3.connect(db_file)  
cur = conn.cursor()    
cur.execute("drop table if exists sample_pop_link")
cur.execute("""create table sample_pop_link
               (spl_id int primary key,
                pop varchar(50),
                sample_id varchar(50))""")

cur.execute(f"create index idx_spl_pop on sample_pop_link(pop)")
cur.execute(f"create index idx_spl_sample_id on sample_pop_link(sample_id)")


#process popfile and identify the samples associated with the pop
with open(pop_file, 'r') as s:
  for line in s:
    line = line.strip()
    pop, species, locations = line.split("\t")
    
    #pop defined for all locations for that species if no location specified
    if locations == 'ALL':
      cur.execute(f"""insert into sample_pop_link
                      select null, '{pop}', sample_id
                      from sample_species
                      where species = '{species}'""")      
    else:
      location_str = "'" + "', '".join(locations.split('~')) + "'"
      
      cur.execute(f"""insert into sample_pop_link
                      select null, '{pop}', sample_id
                      from sample_species
                      where species = '{species}'
                      and location in ({location_str})""")

conn.commit()
conn.close()        