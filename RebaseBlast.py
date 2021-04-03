#import os to work with UNIX commands
import os

import urllib
from Bio import Entrez
import gzip
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-e', '--email', help="Enter your email for Entrez here", required=True)
parser.add_argument('-i', '--input', help="Input file name (list of NCBI accession numbers, one on each line)", required=True)
args = parser.parse_args()

Entrez.email = args.email

#import csv and operator for working with csv files
import csv
import operator

from Bio import SeqIO

#local blast DB is named LocalRebaseDB

sequences='formattedSeqs.fasta'
print("Reading REBASE sequences from " + sequences + "...")
records=list(SeqIO.parse(sequences,'fasta'))

output_name = "output.txt"
final_output=open(output_name,'w')
fileCount = 0

ids=[]
for line in open(args.input).readlines():
  ids.append(line.strip())

os.system('mkdir sequences')

idsused = 1
for id in ids:
  print("Downloading " + id + " from NCBI... (" + str(idsused) + ")")
  handle = Entrez.esearch(db="assembly", term=id+"[id]", retmax='1')
  record = Entrez.read(handle)
  searchid=record['IdList'][0]
  esummary_handle=Entrez.esummary(db="assembly",id=searchid, report="full")
  esummary_record = Entrez.read(esummary_handle)
  url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
  label = os.path.basename(url)
  link = os.path.join(url,label+'_genomic.fna.gz')
  link=link.replace('\\','/')
  urllib.request.urlretrieve(link, 'sequences/' + f'{label}.fna.gz')
  idsused+=1

for fileName in os.listdir("sequences"):

  with gzip.open("sequences/" + fileName,'rb') as f_in:
    with open("sequences/" + fileName[:fileName.index('.gz')],'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)

  temp_output_file='temp.csv'
  
  fileCount += 1
  
  print("Running " + fileName + "..." + "(" + str(fileCount) + ")")
  
  #This blast command will generate a csv formatted output file containing the Query Seq-id, the subject seq-id, the e-value, the query coverage, the percent identity, the bitscore, and the alignment length for each hit.
  blast_command = 'blastn -ungapped -query sequences/' + fileName[:fileName.index('.gz')] + ' -db LocalRebaseDB -out ' + temp_output_file + ' -outfmt "10 qseqid sseqid evalue qcovs pident bitscore length"'
  os.system(blast_command)
  
  #read the csv file and sort it by bitscore
  with open(temp_output_file) as csvFile:
    reader = csv.reader(csvFile, delimiter=',')
    sortedList = sorted(reader,key=lambda x: float(x[5]),reverse=True) #blast output is sorted by bitscore 
  
  topHit=sortedList[0] #now we have just the top hit for the multi-fasta assembly
  
  topHitName=topHit[1]
  for record in records:
    if topHitName == record.id:
      description = record.description
      topHitSequence = str(record.seq)
      
  hitInfo = description.split("\t")
  hitDict = {}
  for info in hitInfo:
    hitDict[info[:info.index(":")]] = info[info.index(":")+1:]
  
  genomeName=fileName #genome name
  contig=topHit[0] #contig
  nameOfSystem = hitDict["REBASE"]
  systemType = hitDict["EnzType"]
  
  final_output.write("Genome file: " + genomeName + "\n")
  final_output.write("Contig: " + contig + "\n")
  final_output.write("Name of System: " + nameOfSystem + "\n")
  final_output.write("System Type: " + systemType + "\n")
  final_output.write(topHitSequence + "\n\n")

print("Cleaning up...")
os.system("rm temp.csv")
print("Done! Output can be found in " + output_name + ".")
#output needed: Genome name - assembly file name
#contig - provided by qseqid in csv file
#name of system - should be provided in top hit name
#system type - should be provided in top hit name
#closest sequence- can get from formatted seqs file