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
  * os
  * csv
  * operator
  * sys

# Instructions

## Using Test Data:
