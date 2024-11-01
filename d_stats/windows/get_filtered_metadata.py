import sys
import gzip
import sqlite3

meta_in = str(sys.argv[1])
vcf_file = str(sys.argv[2])
db_file = str(sys.argv[3])

#create sample_species table
conn = sqlite3.connect(db_file)  
cur = conn.cursor()    
cur.execute("drop table if exists sample_species")
cur.execute("""create table sample_species
               (sample_id varchar(100) primary key,
                species varchar(50),
                location varchar(50),
                vcf_order int)""")


#define short species names which will also serve as pop names
short_species = {'melanogaster': 'mel', 
              'simulans': 'sim',
              'sechellia': 'sech',
              'sim_sech_hybrid': 'ssh',
              'mauritiana': 'mau'}


#process existing metadata file
species = {'D_MELANOGASTER': 'mel'}
locations = {'D_MELANOGASTER': 'Unknown'}
with open(meta_in, 'r') as m:
  next(m)
  for line in m:
    line = line.strip()
    sample, spec, location, date = line.split("\t")
    species[sample] = short_species[spec]
    if location == '':
      location = 'Unknown'
    locations[sample] = location


#insert filtered metadata into sample_species
cur_insert = conn.cursor()    
vcf_order = 0
with gzip.open(vcf_file, 'rt') as v: 
  for line in v:
    if line[0] == '#' and line[:2] != '##':
      line = line.strip()
      values = line.split('\t')
      for sample in values[9:]:
        cur_insert.execute(f"""insert into sample_species
                               values
                               ('{sample}', '{species[sample]}', '{locations[sample]}', {vcf_order})""")
        vcf_order += 1
      break
      
conn.commit()
conn.close()        