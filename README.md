# RM Systems BLAST
COMP 383/483 Project

Authors: Thomas Scott, Anne Jankowski, and Matthew Loffredo

This tool is designed to build upon the already existing BLAST tool for REBASE (the Restriction Enzyme Database), located at http://tools.neb.com/blast/, allowing sequences of any length and batch runs of multiple genomes. The output provided (CSV and FASTA) can then be easily used for additional analysis or exploration. It consists of two python scripts, RebaseSetup and RebaseBlast. RebaseSetup retrieves sequence files from REBASE's FTP site, formats them, and creates a local BLAST+ database. These files are updated monthly (see http://rebase.neb.com/rebase/rebase.files.html), so this script only needs to be run when the most updated files are needed. RebaseBlast can then BLAST multiple genomes/assemblies against this database (ungapped), and return information about each hit to two output files: CSV output containing BLAST results, and a FASTA output with information from REBASE.

# Requirements

_This tool was made in and meant to be run in a Linux or OSX terminal_

* Blast+ command line application (https://www.ncbi.nlm.nih.gov/books/NBK279690/)
* Python3 and packages used:
  * Biopython (https://biopython.org/wiki/Download)
    * SeqIO
    * Entrez
  * os
  * csv
  * operator
  * sys
  * urllib
  * FTP
  * gzip
  * shutil
  * argparse

# Instructions

### RebaseSetup.py

In order to run the tool, the RebaseSetup script must be run first. This downloads required files and sets up the BLAST database. RebaseSetup only needs to be run once initially. If you have not used the tool in a while, it is recommended to run it again as REBASE files may have been updated.

**Command for running RebaseSetup:**
```
python3 RebaseSetup.py --email EMAIL
```
Where EMAIL is your email (for accessing REBASE's FTP site). This is a required parameter.

Again, this script only needs to be run once initially, and only when REBASE files are updated afterwards.

### RebaseBlast.py

Once you have run RebaseSetup, you are ready to BLAST sequences against REBASE. You can provide input in a couple ways. First, you could provide a text file with NCBI genome accession numbers (i.e. GCF_001641215.1) listed one per line (see AccessionNumbers10.txt for an example file) using the parameter --input. The script will download the assembly FASTA files for you to a directory called 'sequences' and unzip them. Alternatively, you may manually place any already downloaded FASTA files into a 'sequences' directory yourself. You may use both manually downloaded sequences and accession numbers as well, if you so choose. In that case, the script will be sure not to download duplicates if they are already in the 'sequences' directory.

**Command for running RebaseBlast:**
```
python3 RebaseBlast.py --email EMAIL --input INPUT
```

Where EMAIL is your email (used for Entrez), and INPUT is the input file explained above. Your email is a required parameter.

## Output

Two output files will be produced: one FASTA formatted (RMBlastFasta.fasta) and one CSV formatted (RMBlastCSV.csv). The CSV file contains the top hits with BLAST information for all the genomes ran (based on bitscore), and is good for getting a general picture of the BLAST results. If you would like to know more about the RM systems matched to your input, their sequences and protein IDs are found in the FASTA file. Furthermore, full BLAST results for each genome ran are provided in the directory 'full_blast_results'. These outputs are helpful on their own, but we encourage you to use them for further analysis! An example of analysis you could do is provided, and is explained below.

## Using Test Data:

Test Data (10 accession numbers) is provided in AccessionNumbers10.txt. To run this, run these commands:

```
python3 RebaseSetup.py --email EMAIL
```

```
python3 RebaseBlast.py --email EMAIL --input AccessionNumbers10.txt
```

### Further Analysis Example

The output from this tool is great for performing additional analysis. For example, you may compare the presence of RM systems in your genomes to the presence of plasmids or phages. In this example, plasmid and phage presence was analyzed in the genomes from AccessionNumbers10.txt and AccessionNumbersAll.txt, with the results for those being found in plasmid_finder_results_putonti.tsv and phage_results_putonti.tsv. Two additional python scripts (ParsePlasmidGenomes and ParsePhageGenomes) were used to compare these results to the RM systems BLAST results, giving a picture of what bacteria with RM system types may also have plasmids and phages. It also goes further to try and answer the question: for each RM system type, which genera that contain them are more likely to have plasmids or phages? You may run these scripts with the test data to see this example.

To do this, use these two commands:

```
python3 ParsePlasmidGenomes.py
```

```
python3 ParsePhageGenomes.py
```

Output for these scripts will be found in plasmidCompareResults.csv and phageCompareResults.csv.
