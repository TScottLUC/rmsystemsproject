#import os to work with UNIX commands
import os

#assume the local blast db is already created - for example with title 'rebase'
#for now, this code assumes there is only one multi-FASTA file being used
#it will later need to BLAST many multi-FASTA files in one run

input_file = "example_file.fasta"
output_file = "output.csv"

#This blast command will generate a csv formatted output file containing the Query Seq-id, the subject seq-id, the e-value, the query coverage, the percent identity, the score, and the alignment length for each hit.
blast_command = 'blastn -query ' + input_file + ' -db rebase -out ' + output_file + ' -outfmt 10 qseqid sseqid evalue qcovs pident score length'
os.system(blast_command)

#from this csv file, we will need to parse just the top hit for each multi-fasta assembly run

#output needed: genome name (which genome the hit came from) - can get this from input file
#contig - may be provided from qseqid
#name of system - should be provided in top hit name
#system type - should be provided in top hit name
#closest sequence- can get from the dna_seqs file?