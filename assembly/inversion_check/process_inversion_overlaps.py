import sqlite3
import os
import tempfile
import sys
import re
import glob
import argparse

#python3 process_inversion_overlaps.py inversion_overlaps_contigs/*inv_overlap_contigs.txt -l contigs trimmedReads


def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("seq_types", metavar='N', type=str, nargs='*')

  
  return parser.parse_args()


def main():
  args = parse_args()
  seq_types = args.seq_types  
  
  
  db_file = 'minimap_results.db'
  conn = sqlite3.connect(db_file)
    
  cur = conn.cursor()                            

  cur.execute("drop table if exists minimap_results")
  cur.execute("""create table minimap_results
                 (mr_id int primary key,
                  sample_id varcahr(50),
                  inv_name varchar(20),
                  seq_type varchar(12),
                  contig_name varchar(50),
                  contig_len int,
                  contig_start int,
                  contig_end int,
                  chrom varchar(5),
                  start int,
                  end int,
                  strand varchar(1))""")
   
  cur.execute("create index idx_mr_sample_id on minimap_results(sample_id)")
  cur.execute("create index idx_mr_inv_name on minimap_results(inv_name)")
  cur.execute("create index idx_mr_read_name on minimap_results(contig_name)")
  cur.execute("create index idx_mr_start on minimap_results(start)")
  conn.close()

  
  mr_id = 0
  for seq_type in seq_types:
    with tempfile.NamedTemporaryFile(mode='w') as t:
      for inv_overlap_file in glob.glob(os.path.join('inversion_overlaps_' + seq_type, '*_inv_overlap_' + seq_type + '.txt')):
        
        m = re.search(r'.+/(.+)_inv_overlap.+.txt', inv_overlap_file)
        sample_id = m.group(1)
        
        with open(inv_overlap_file, 'r') as i: 
          for line in i:
            mr_id += 1
            line = line.strip()
            values = line.split('\t')
                    
            t.write(f"{mr_id}\t{sample_id}\t{values[0]}\t{seq_type}\t{values[1]}\t{values[2]}\t{values[3]}\t{values[4]}\t{values[5]}\t{values[6]}\t{values[7]}\t{values[8]}\n")   
                             
      t.flush()      
      os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} minimap_results" """)


if __name__ == '__main__':
  main()
