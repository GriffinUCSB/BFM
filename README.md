# BFM
A collection of python scripts that help with bioinformatics file management
all scripts require manual editing to change filepaths and output names, this will be edited in future versions

## create_custom_db.py
creates a custom library file given a count table and a list of library files. Outputs .txt with columns sgID, sublibrary, gene, transcripts, sequence.

## create_wide_ct.py
combines a list of count tables to be used for MAGeCK.

#### Usage

*python create_wide_ct.py [Options] -n [file names]*

##### File names
input file names in specific order: output_file.txt count_table1.txt ... count_tableN.txt

##### Options 
*--header*    
this option should be used when creating a new output file. Adds header. Do not use if appending to previously created output.   

## find_unmapped.py
take .txt of unmapped reads and match them to a list of libraries. Unmapped reads file should be created with the terminal command below. 
samtools view -F 4 file.bam | awk '{print $10}' | sort | uniq -c | sort -k1nr | head -n 20 > unmapped.txt
Can be used to increase alignment rates by adding unmapped reads to .fa and recreating alignment lib.

#### Usage

*python find_unmapped.py [Options] -n [file names]*

##### File names
input file names in specific order:
unmapped_reads.txt output_file.txt library1.txt library2.txt ... libraryN.txt

##### Options
*--header*     
this option should be used when creating a new output file. Adds header. Do not use if appending to previously created output.       
*-p, --percentage*          
provides a percent complete marker in terminal.    
*-i, --include-unfound*     
includes reads that remain unfound even after the search is completed.   

## merge_gene_tables.py
takes two genetables in .xls format (from screen_processing script) and merges them. When configuring the longer table should go first. Be sure the longer table contains all of the genes present in the shorter table for script to work. 

## proteomics_to_sgid
takes a proteomics xlsx and turns it into a list of single guide RNAs

