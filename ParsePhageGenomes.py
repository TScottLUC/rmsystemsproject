# parse accession numbers of genomes with phages into txt file
import csv

fname = "phage_results_putonti.tsv"

#parse the info from phage results to an easier to use file
print("Parsing genome accessions with phages to file phageAccessions.txt...")
with open(fname) as f_in, open('phageAccessions.txt', 'w') as f_out:
  rd = csv.reader(f_in, delimiter="\t", quotechar='"')
  headers = next(rd)
  for row in rd:
    f_out.write(row[0]+','+row[1]+','+row[2]+'\n')
  f_in.close(); f_out.close();
  
  
phageResultsFileName = 'phageCompareResults.csv'
print("Running phage analysis... (output will be in " + phageResultsFileName + ")")

#read the phage results, putting the information into a dictionary for each accession
accessionDict = {}
for line in open('phageAccessions.txt'):
  accessionInfo = line.split(',')
  accessionDict[accessionInfo[0]] = []
  accessionDict[accessionInfo[0]].append(accessionInfo[1].split(" ")[0])
  accessionDict[accessionInfo[0]].append(accessionInfo[2].strip()) 
  #key = accession, value = [genus, phage #], RM system type will be appended later

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
  accessionDict[accession].append(enzType) #RM system type is now in the information for each accession (also contains genus and # for phage)  

#open an output file for the phage results
phageResults = open(phageResultsFileName,'w')

#initialize two dictionaries: one to holds the number of 'yes' for type and one that holds all for type
typeYes = {}
typeAll = {}

#go through each accession that was ran in the BLAST analysis, counting the number of each
#RM system type, and the number of that type that have phages as well
for accession in accessionsWithBlastInfo:
  rmtype = accessionDict[accession][2]
  hasPhage = (int(accessionDict[accession][1]) > 0)
  
  #add to total of RM System type
  if rmtype not in typeAll:
    typeAll[rmtype] = 1
  else:
    typeAll[rmtype] += 1
    
  #add to number that has phages if applicable
  if hasPhage:
    if rmtype not in typeYes:
      typeYes[rmtype] = 1
    else:
      typeYes[rmtype] += 1

#write headers for the type information to be written
phageResults.write('RM System Type,Total # of System Type,# with Phage(s),% with Phage(s)\n')
for rmType in typeAll.keys():
  
  phageResults.write(rmType + ',') #write the type name
  phageResults.write(str(typeAll[rmType]) + ',') #write the total number of that type
 
  #write the number of 'yes' values for that type
  try:
    phageResults.write(str(typeYes[rmType]) + ',')
  except KeyError: #if it has no 'yes' values, write 0
    typeYes[rmType] = 0
    phageResults.write(str(typeYes[rmType]) + ',')
    
  #write the % of type that have phage
  phageResults.write(str(float((typeYes[rmType]/typeAll[rmType]))*100) + '\n')
phageResults.write('\n')

#for each type (already written to the file), analyze which genera containing that type also have phages
#are certain genera more likely to have or not have phages?
for type in typeYes:

  #initialize dictionaries to hold the number of 'yes' for each genera and the total for each genera
  genusYes = {}
  genusTotal = {}
  
  #write the type and headers to the file
  phageResults.write(type + '\n')
  phageResults.write('Genus,Total # of Genus,# with Phage,% with Phage\n')
  
  #go through the accessions
  for acc in accessionsWithBlastInfo:
  
    accInfo = accessionDict[acc]
    
    #extract the RM type that was found from that accession
    #and check if it matches the current type
    rmtype = accInfo[2]
    if rmtype == type:
      
      #if it matches the current type, extract the genus and check if it has a phage
      genus = accInfo[0]
      hasPhage = (int(accInfo[1]) > 0)
      
      #add the genus to the total
      if genus not in genusTotal:
        genusTotal[genus] = 1
      else:
        genusTotal[genus] += 1
        
      #add the genus to the number with phages, if applicable
      if hasPhage:
        if genus not in genusYes:
          genusYes[genus] = 1
        else:
          genusYes[genus] += 1
  
  #for each genus under the current type, write its info to the file
  for genus in genusTotal.keys():
  
    #write the genus
    phageResults.write(genus + ',')
    
    #write the total number of that genus with the RM system type
    phageResults.write(str(genusTotal[genus]) + ',')
    
    #write the total number of that genus with phage
    try:
      phageResults.write(str(genusYes[genus]) + ',')
    except KeyError:
      genusYes[genus] = 0
      phageResults.write(str(genusYes[genus]) + ',')
    
    #write the % of genus with phage
    phageResults.write(str(float((genusYes[genus]/genusTotal[genus]))*100) + '\n')
  phageResults.write('\n')