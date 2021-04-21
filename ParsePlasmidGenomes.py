# parse accession numbers of genomes with plasmids into txt file
import csv

fname = "plasmid_finder_results_putonti.tsv"

#parse the info from phage results to an easier to use file
print("Parsing genome accessions with plasmids to file plasmidAccessions.txt...")
with open(fname) as f_in, open('plasmidAccessions.txt', 'w') as f_out:
  rd = csv.reader(f_in, delimiter="\t", quotechar='"')
  headers = next(rd)
  for row in rd:
    f_out.write(row[1]+','+row[2]+','+row[3]+'\n')
  f_in.close(); f_out.close();


plasmidResultsFileName = 'plasmidCompareResults.csv'
print("Running plasmid analysis... (output will be in " + plasmidResultsFileName + ")")

#read the plasmid results, putting the information into a dictionary for each accession
accessionDict = {}
for line in open('plasmidAccessions.txt'):
  accessionInfo = line.split(',')
  accessionDict[accessionInfo[0]] = []
  accessionDict[accessionInfo[0]].append(accessionInfo[1].split(" ")[0])
  accessionDict[accessionInfo[0]].append(accessionInfo[2].strip()) 
  #key = accession, value = [genus, plasmid yes/no], RM system type will be appended later
  
#initialize a list for accessions with BLAST info: if it doesn't have BLAST info, it won't be used in the analysis
accessionsWithBlastInfo = []
#read all lines of the CSV (blast info) file, skipping the first line
for line in open('RMBlastCSV.csv').readlines()[1:]:
  
  #extract the accession number from the file name
  blastInfo = line.split(',')
  accession = blastInfo[0]
  accession = accession[:accession.index('_ASM')]
  accessionsWithBlastInfo.append(accession)
    
  #extract the enzyme type that was found
  enzType = blastInfo[3]
  accessionDict[accession].append(enzType) #RM system type is now in the information for each accession (also contains genus and yes/no for plasmid)  
  
#open an output file for the plasmid results
plasmidResults = open(plasmidResultsFileName,'w')

#initialize two dictionaries: one to holds the number of 'yes' for type and one that holds all for type
typeYes = {}
typeAll = {}

#go through each accession that was ran in the BLAST analysis, counting the number of each
#RM system type, and the number of that type that have plasmids as well
for accession in accessionsWithBlastInfo:
  rmtype = accessionDict[accession][2]
  hasPlasmid = (accessionDict[accession][1] == 'yes')
  
  #add to total of RM System type
  if rmtype not in typeAll:
    typeAll[rmtype] = 1
  else:
    typeAll[rmtype] += 1
    
  #add to number that has plasmids if applicable
  if hasPlasmid:
    if rmtype not in typeYes:
      typeYes[rmtype] = 1
    else:
      typeYes[rmtype] += 1

#write headers for the type information to be written
plasmidResults.write('RM System Type,Total # of System Type,# with Plasmid,% with Plasmid\n')
for rmType in typeAll.keys():

  plasmidResults.write(rmType + ',') #write the type name
  plasmidResults.write(str(typeAll[rmType]) + ',') #write the total number of that type
  
  #write the number of 'yes' values for that type
  try:
    plasmidResults.write(str(typeYes[rmType]) + ',')
  except KeyError: #if it has no 'yes' values, write 0
    typeYes[rmType] = 0
    plasmidResults.write(str(typeYes[rmType]) + ',')
    
  #write the % of type that have plasmid
  plasmidResults.write(str(float((typeYes[rmType]/typeAll[rmType]))*100) + '\n')
plasmidResults.write('\n')

#for each type (already written to the file), analyze which genera containing that type also have plasmids
#are certain genera more likely to have or not have plasmids?
for type in typeYes:

  #initialize dictionaries to hold the number of 'yes' for each genera and the total for each genera
  genusYes = {}
  genusTotal = {}
  
  #write the type and headers to the file
  plasmidResults.write(type + '\n')
  plasmidResults.write('Genus,Total # of Genus,# with Plasmid\n')
  
  #go through the accessions
  for acc in accessionsWithBlastInfo:
  
    accInfo = accessionDict[acc]
    
    #extract the RM type that was found from that accession
    #and check if it matches the current type
    rmtype = accInfo[2]
    if rmtype == type:
    
      #if it matches the current type, extract the genus and check if it has a plasmid
      genus = accInfo[0]
      hasPlasmid = (accInfo[1] == 'yes')
      
      #add the genus to the total
      if genus not in genusTotal:
        genusTotal[genus] = 1
      else:
        genusTotal[genus] += 1
        
      #add to the number with plasmids, if applicable
      if hasPlasmid:
        if genus not in genusYes:
          genusYes[genus] = 1
        else:
          genusYes[genus] += 1
          
  #for each genus under the current type, write its info to the file
  for genus in genusTotal.keys():
    
    #write the genus
    plasmidResults.write(genus + ',')
    
    #write the total number of that genus with the RM system type
    plasmidResults.write(str(genusTotal[genus]) + ',')
    
    #write the total number of that genus with plasmid
    try:
      plasmidResults.write(str(genusYes[genus]) + ',')
    except KeyError:
      genusYes[genus] = 0
      plasmidResults.write(str(genusYes[genus]) + ',')
      
    #write the % of genus with plasmid
    plasmidResults.write(str(float((genusYes[genus]/genusTotal[genus]))*100) + '\n')
  plasmidResults.write('\n')    