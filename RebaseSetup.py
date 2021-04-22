#import os for working with unix commands
import os

#import SeqIO for parsing FASTA files
from Bio import SeqIO

#import FTP for using REBASE's FTP site
from ftplib import FTP

#import argparse for command-line parameters
import argparse

#create the command-line parameter for email (for REBASE FTP)
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--email', help="Enter your email for REBASE FTP site here", required=True)
args = parser.parse_args()
email = args.email

#retrieve the REBASE DNA sequence files using FTP
print('Retrieving sequences from REBASE')
ftp = FTP('ftp.neb.com')
ftp.login(passwd=args.email)
ftp.cwd('pub/rebase')
# ftp.retrlines('LIST')
with open('dna_seqs.txt', 'wb') as seqs_file:
  ftp.retrbinary('RETR dna_seqs.txt', seqs_file.write)
  seqs_file.close()

# clean up lines from file received from REBASE
print('Formatting sequences received from REBASE')
with open('dna_seqs.txt') as infile, open('rebase_seqs.fasta', 'w') as outfile:
  for line in infile:
    if not line.strip(): continue  # skip the empty line
    #if 'ThisgenbanknumberdoesnothaveaDNAsequenceassociatedwithit.' in line: continue # skip line w/ no seq 
    line = line.replace("<>", "") # remove <> from lines
    outfile.write(line)  # non-empty line. Write it to output
  infile.close(); outfile.close()

# writes to file in fasta format, makeblast command fails without this
with open('rebase_seqs.fasta', "r") as input_handle, open("formattedSeqs.fasta", "w") as output_handle:
  sequences = SeqIO.parse(input_handle, "fasta")
  SeqIO.write(sequences, output_handle, "fasta")
  output_handle.close()
  input_handle.close()
  
# Delete temp. seqs file
os.system("rm -f rebase_seqs.fasta")

# Make local blast DB
print('Making local BLAST db')
input_file = 'formattedSeqs.fasta'
makeblast_command = "makeblastdb -in "+input_file+" -out LocalRebaseDB -title LocalRebaseDB -dbtype nucl"
os.system(makeblast_command)

