# parse accession numbers of genomes with plasmids into txt file
import csv

fname = "phage_results_putonti.tsv"

print("Parsing genome accessions with phages to file phageAccessions.txt...")
with open(fname) as f_in, open('phageAccessions.txt', 'w') as f_out:
  rd = csv.reader(f_in, delimiter="\t", quotechar='"')
  headers = next(rd)
  for row in rd:
    if(int(row[2]) > 0):
      f_out.write(row[0]+', '+row[1]+'\n')
  f_in.close(); f_out.close();
