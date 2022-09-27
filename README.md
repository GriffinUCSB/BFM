# BFM
A collection of python scripts that help with bioinformatics file management
all scripts require manual editing to change filepaths and output names, this will be edited in future versions

## create_custom_db.py
creates a custom library file given a count table and a list of library files. Outputs .txt with columns sgID, sublibrary, gene, transcripts, sequence.

## create_wide_ct.py
combines a list of count tables to be used for MAGeCK.

## find_unmapped.py
take .txt of unmapped reads and match them to a list of libraries. Unmapped reads file should be created with the terminal command below. 
samtools view -F 4 file.bam | awk '{print $10}' | sort | uniq -c | sort -k1nr | head -n 20 > unmapped.txt
Can be used to increase alignment rates by adding unmapped reads to .fa and recreating alignment lib.

## merge_gene_tables.py
takes two genetables in .xls format (from screen_processing script) and merges them. When configuring the longer table should go first. Be sure the longer table contains all of the genes present in the shorter table for script to work. 

