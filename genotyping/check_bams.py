import sqlite3
import os


#python3 check_bams.py




def main():
  db_file = '/work/users/d/t/dturissi/drosophila/ssh/introgression/d_stats/windows/ssh_d_win.db'
  bam_dir = '/proj/matutelb/users/stuckert/simulans/bams'
  
  
  conn = sqlite3.connect(db_file)
      
  sample_query = conn.execute(f"""select sample_id, species
                                  from sample_species""")
  

  for sample_id, species in sample_query:
    bam_file = os.path.join(bam_dir, sample_id + '.dedupd.bam')
    
    if not os.path.exists(bam_file):
      print(sample_id, species)



if __name__ == '__main__':
  main()
