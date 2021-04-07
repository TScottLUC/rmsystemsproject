# parse accession numbers of genomes with plasmids into txt file
import csv

fname = "plasmid_finder_results_putonti.tsv"

print("Parsing genome accessions with plasmids to file plasmidAccessions.txt...")
with open(fname) as f_in, open('plasmidAccessions.txt', 'w') as f_out:
  rd = csv.reader(f_in, delimiter="\t", quotechar='"')
  for row in rd:
    if(row[3] == "yes"):
      f_out.write(row[0]+'\n')
  f_in.close(); f_out.close();
