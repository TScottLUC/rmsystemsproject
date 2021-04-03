# RM Systems BLAST
COMP 383/483 Project

Authors: Thomas Scott, Anne Jankowski, and Matthew Loffredo

This tool is designed to build upon the already existing BLAST tool for REBASE (the Restriction Enzyme Database), located at http://tools.neb.com/blast/, allowing sequences of any length, multi-FASTA file inputs, and batch runs. It consists of two python scripts, RebaseSetup and RebaseBlast. RebaseSetup retrieves sequence files from REBASE's FTP site, and creates a local BLAST+ database. These files are updated monthly (see http://rebase.neb.com/rebase/rebase.files.html), so this script only needs to be run when the most updated files are needed. RebaseBlast can then BLAST multiple genomes/assemblies against this database, and return information about each hit to an output file, including the contig of the hit, the name of the RM system, the type of system, and the sequence from REBASE.

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
  * gzip
  * shutil
  * argparse

# Instructions

Run both Python scripts (RebaseSetup and RebaseBlast), one after the other. RebaseSetup only needs to be run once initially. If you have not used the tool in a while, it is recommended to run it again as REBASE files may have been updated.

### RebaseSetup.py

**Command for running RebaseSetup:**
```
python3 RebaseSetup.py
```

Again, this script only needs to be run once initially, and only when REBASE files are updated afterwards.

### RebaseBlast.py

Input for this tool is a text file with NCBI genome accession numbers (i.e. GCF_001641215.1) listed one per line (see AccessionNumbers10.txt for an example file).

**Command for running RebaseBlast:**
```
python3 RebaseBlast.py --email EMAIL --input INPUT
```

Where EMAIL is your email (used for Entrez), and INPUT is the input file explained above. These are required parameters.

Output will be found in a file named output.txt.

## Using Test Data:

Test Data (10 accession numbers) is provided in AccessionNumbers10.txt. To run this, simply use this command:

```
python3 RebaseBlast.py --email EMAIL --input AccessionNumbers10.txt
```
