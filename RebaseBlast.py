#import os to work with UNIX commands
import os

#import SeqIO for parsing fasta files
from Bio import SeqIO

#import Entrez and urllib to download assembly files from NCBI's FTP
from Bio import Entrez
import urllib

#import gzip and shutil for working with zipped assembly files
import gzip
import shutil

#import argparse for command-line parameters
import argparse

#create the command-line parameters for email (for Entrez) and input file
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--email', help="Enter your email for Entrez here", required=True)
parser.add_argument('-i', '--input', help="Input file name (list of NCBI accession numbers, one on each line)", required=True)
parser.add_argument('-p', '--plasmids', help="Use this tag if you would also like to compare the BLAST results to Plasmid Finder results (from plasmid_finder_results_putonti.tsv)", action="store_true")
args = parser.parse_args()

#set Entrez's email to the email provided
Entrez.email = args.email

#import csv and operator for working with csv files
import csv
import operator

#read REBASE sequence files, as these records will need to be parsed for
#information on each top hit.
sequences='formattedSeqs.fasta'
print("Reading REBASE sequences from " + sequences + "...")
rebase_records=list(SeqIO.parse(sequences,'fasta'))

#read the genome accession #'s from the input file
#(specified in args.input parameter)
ids=[]
for line in open(args.input).readlines():
  ids.append(line.strip())

#make a sequences directory to store all sequence fasta files in
if (os.path.isdir('sequences')):
  os.system('rm -r sequences')
os.system('mkdir sequences')

#this for loop takes each id and downloads it's assembly file from NCBI through FTP
idsused = 0 #count how many files have been downloaded so the user can follow the progress of the script
for id in ids:
  
  idsused+=1
  print("Downloading " + id + " from NCBI... (" + str(idsused) + ")")
  
  #search NCBI assembly database by genome accession # (id), and retrieve only 1 record
  handle = Entrez.esearch(db="assembly", term=id+"[id]", retmax='1')
  record = Entrez.read(handle)
  
  #get the assembly search id for the genome accession provided
  searchid=record['IdList'][0]
  
  #get the assembly summary/report using the assembly search id
  esummary_handle=Entrez.esummary(db="assembly",id=searchid, report="full")
  esummary_record = Entrez.read(esummary_handle,validate=False)
  
  #get the FTP link for the assembly file
  try:
    url = esummary_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
  except IndexError:
    continue
  
  #download the assembly file via FTP to the sequences directory
  label = os.path.basename(url)
  link = os.path.join(url,label+'_genomic.fna.gz')
  link=link.replace('\\','/') #reformat some slashes from the link
  urllib.request.urlretrieve(link, 'sequences/' + f'{label}.fna.gz')

#open the output files
csv_output_name = "RMBlastCSV.csv"
csv_output=open(csv_output_name,'w')
fasta_output_name="RMBlastFasta.fasta"
fasta_output=open(fasta_output_name,'w')

#write the headers to the simple output file
csv_output.write('Genome,Contig,RMSystemName,SystemType,E-Value,% Identity,Query Coverage,Bitscore,Length\n')

#count the files BLASTed so the user can follow the progress of the script
fileCount = 0

#this for loop goes through each file in the sequences directory, unzips it, and BLASTs it against the local REBASE database created in the RebaseSetup script. It then writes information about the top hit to the output file.
for fileName in os.listdir("sequences"):

  unzippedFileName = fileName[:fileName.index('.gz')]
  #unzip the .gz assembly files
  with gzip.open("sequences/" + fileName,'rb') as f_in:
    with open("sequences/" + unzippedFileName ,'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)
      
  #delete the original zipped file
  os.system('rm sequences/' + fileName)

  #temporarily holds the BLAST results for each file
  temp_output_file='temp.csv'
  
  fileCount += 1
  print("Running " + unzippedFileName + "..." + "(" + str(fileCount) + ")")
  
  #This blast command will generate a csv formatted output file containing the Query Seq-id, the subject seq-id, the e-value, the query coverage, the percent identity, the bitscore, and the alignment length for each hit.
  #the local BLAST db created in RebaseSetup is named LocalRebaseDB
  blast_command = 'blastn -ungapped -query sequences/' + fileName[:fileName.index('.gz')] + ' -db LocalRebaseDB -out ' + temp_output_file + ' -outfmt "10 qseqid sseqid evalue pident qcovs bitscore length"'
  os.system(blast_command)
  
  #read the csv file and sort it by bitscore
  with open(temp_output_file) as csvFile:
    reader = csv.reader(csvFile, delimiter=',')
    sortedList = sorted(reader,key=lambda x: float(x[5]),reverse=True) #blast output is sorted by bitscore 
  
  topHit=sortedList[0] #now we have just the top hit for the multi-fasta assembly
  
  #sseqid is in index 1
  topHitName=topHit[1]
  
  #find the top hit in REBASE's records
  for record in rebase_records:
    if topHitName == record.id:
      description = record.description #the description contains several pieces of information about the RM system
      topHitSequence = str(record.seq) #get the sequence of the top hit
      
  #the description is tab-delimitted, so split it to get each portion
  hitInfo = description.split("\t")
 
  #create a dictionary to store each piece of information
  #(the description follows the format KEY:VALUE\tKEY:VALUE\t, etc.)
  hitDict = {}
  for info in hitInfo:
    hitDict[info[:info.index(":")]] = info[info.index(":")+1:]
  
  genomeName=unzippedFileName #genome name is the file name
  contig=topHit[0] #contig where the top hit matched to is provided in index 0 of the temporary csv output
  
  #get the name of the system and its type from the REBASE description
  nameOfSystem = hitDict["REBASE"]
  systemType = hitDict["EnzType"]
  try:
    proteinID = hitDict["ProteinId"]
  except KeyError:
    proteinID = "N/A"
  
  eValue = topHit[2]
  pIdent = topHit[3]
  qCov = topHit[4]
  bitscore = topHit[5]
  length = topHit[6]
  
  #write relevant output to the output files
  csv_output.write(genomeName + ",")
  csv_output.write(contig + ",")
  csv_output.write(nameOfSystem + ",")
  csv_output.write(systemType + ",")
  csv_output.write(eValue + ",")
  csv_output.write(pIdent + ",")
  csv_output.write(qCov + ",")
  csv_output.write(bitscore + ",")
  csv_output.write(length + "\n")
  
  fasta_output.write(">Genome: " + genomeName + "\t")
  fasta_output.write("Contig: " + contig + "\t")
  fasta_output.write("RMSystemName: " + nameOfSystem + "\t")
  fasta_output.write("SystemType: " + systemType + "\t")
  fasta_output.write("ProteinID: " + proteinID + "\t")
  fasta_output.write("E-Val: " + eValue + "\t")
  fasta_output.write("Bitscore: " + bitscore + "\n")
  fasta_output.write(topHitSequence + "\n")
  
#close the output files
csv_output.close()
fasta_output.close()

#perform plasmid analysis if --plasmids was used 
if args.plasmids:

  #run the plasmid parser
  os.system('python3 ParsePlasmidGenomes.py')
  
  #add the accession numbers with plasmids to a list
  plasmidYes = []
  for line in open('plasmidAccessions.txt'):
    plasmidYes.append(line.strip())
    
  #initialize two dictionaries: one to holds the number of 'yes' for type and one that holds all for type
  typeYes = {}
  typeAll = {}
  
  #read all lines of the CSV (blast info) file, skipping the first line
  for line in open('RMBlastCSV.csv').readlines()[1:]:
  
    #extract the accession number from the file name
    blastInfo = line.split(',')
    accession = blastInfo[0]
    accession = accession[:accession.index('_ASM')]
    
    #extract the enzyme time that was found
    enzType = blastInfo[3]
    
    #add 1 to the enzyme type in the typeAll dictionary
    if enzType in typeAll:
      typeAll[enzType] += 1
    else:
      typeAll[enzType] = 1
      
    #add 1 to the enzyme type in the typeYes dictionary if it is present in plasmidYes
    #(there has been shown a plasmid with a RM system of this type)
    if accession in plasmidYes:
      if enzType in typeYes:
        typeYes[enzType] += 1
      else:
        typeYes[enzType] = 1
  
  #open an output file for the plasmid results
  plasmidResultsFileName = 'plasmidCompareResults.csv')
  plasmidResults = open(plasmidResultsFileName,'w')
  
  #write headers for the information to be written
  plasmidResults.write('RM System Type,Total # of System Type,# with Plasmid,% with Plasmid\n')
  for rmType in typeAll.keys():
      plasmidResults.write(rmType + ',') #write the type name
      plasmidResults.write(str(typeAll[rmType]) + ',') #write the total number of that type
      #write the number of 'yes' values for that type
      try:
        plasmidResults.write(str(typeYes[rmType]) + ',')
      except KeyError: #if it has no 'yes' values, write 0
        typeYes[rmType] = 0
        plasmidResults.write(str(typeYes[rmType]) + '.,')
      #write the % of type that have plasmid
      plasmidResults.write(str(float((typeYes[rmType]/typeAll[rmType]))*100) + '\n')


print("Cleaning up...")
os.system("rm temp.csv")
print("Done!")
print("Output with more BLAST information (without sequences, csv formatted) can be found in " + csv_output_name + ".")
print("Output with sequences from REBASE and protein IDs (fasta formatted) can be found in " + fasta_output_name + ".")

if args.plasmids:
  print("Output with comparison to the plasmid finder results can be found in " + plasmidResultsFileName + ".")