#for now, this code assumes there is only one multi-FASTA file being used
#it will later need to BLAST many multi-FASTA files in one run

#import os to work with UNIX commands
import os

#import csv and operator for working with csv files
import csv
import operator

#local blast DB is named LocalRebaseDB

input_file='shortened.fasta' #EXAMPLE FILE
output_file='temp.csv'

#This blast command will generate a csv formatted output file containing the Query Seq-id, the subject seq-id, the e-value, the query coverage, the percent identity, the bitscore, and the alignment length for each hit.
blast_command = 'blastn -query ' + input_file + ' -db LocalRebaseDB -out ' + output_file + ' -outfmt "10 qseqid sseqid evalue qcovs pident bitscore length"'
os.system(blast_command)

#read the csv file and sort it by bitscore
with open(output_file) as csvFile:
  reader = csv.reader(csvFile, delimiter=',')
  sortedList = sorted(reader,key=lambda x: float(x[5]),reverse=True) #blast output is sorted by bitscore
  
topHit=sortedList[0] #now we have just the top hit for the multi-fasta assembly

#output needed: genome name (which genome the hit came from) - can get this from input file
#contig - may be provided from qseqid
#name of system - should be provided in top hit name
#system type - should be provided in top hit name
#closest sequence- can get from the dna_seqs file?