import os
import glob
import re


def main():
  base_dir = '/work/users/d/t/dturissi/drosophila/ssh/assembly'
  fastq_dir = os.path.join(base_dir, 'fastqs')
  read_length_file = os.path.join(base_dir, 'sech_fastq_read_lengths.txt')
  
  
  with open(read_length_file, 'w') as o:
    read_lens = {}
    for fastq_file in glob.glob(os.path.join(fastq_dir, '*.fastq')): 
      m = re.search(r'.+/(.+).fastq', fastq_file)
      fly_line = m.group(1)
      print(fly_line)
      
      with open(fastq_file, 'r') as f: 
        for line_1 in f:
          line_2 = next(f)
          line_3 = next(f)
          line_4 = next(f)
          
          read_len = len(line_2.strip())
        
          if read_len in read_lens:
            read_lens[read_len] += 1
          else:
            read_lens[read_len] = 1
        
      for read_len in sorted(read_lens):
        o.write(f"{fly_line}\t{read_len}\t{read_lens[read_len]}\n")


if __name__ == '__main__':
  main()
