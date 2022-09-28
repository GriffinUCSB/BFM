#find unmapped reads V2 6/14/2021
#Griffin Kramer

#!important note! the built in delimiter for unmapped is a space based on the samtools output of the line:
# samtools view -F 4 file.bam | awk '{print $10}' | sort | uniq -c | sort -k1nr | head -n 20 > unmapped.txt
#the delimiter for librarys is a tab ("\t")
#both of these can easily be changed manually below

#imports
import csv, pandas as pd
import os
import argparse
import warnings
import sys

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='proccess user in')
parser.add_argument("--header", help="add header to output file", action="store_true")
parser.add_argument('-p', "--percentage", help="print complete percentsge", action="store_true")
parser.add_argument('-i', "--include_unfound", help="include reads that remain unfound in output", action="store_true")
parser.add_argument('-n', '--file_names',nargs='+', default=[], help='Input.txt file names in this order: unmapped_reads.txt output.txt libfile1.txt libfileN.txt')
args=parser.parse_args()
#variables for configuration
unmapped_file_name = args.file_names[0]
library_file_names = args.file_names[2:]

if args.header:
    header = True #If this program is being used to create a new .csv, set header to True in config
else:
    header = False

def search_libraries(unmapped, library, header):

    #prepare output dataframe
    output_columns = ["gene", "seq", "unmapped_file", "freq", "lib_file"]
    outdf = pd.DataFrame(columns = output_columns)

    #load unmapped sequence .txt file into dataframe
    unmappeddf = pd.read_csv(unmapped, header = None, delimiter="\t")
    unmappeddf.columns = ['freq', 'seq']
    rowUnmapped, colUnmapped = unmappeddf.shape #get shape to use as length for iteration

    #load libraries into dataframes
    librarydfs = []
    for library_iterator in range (len(library)):
        librarydfs.append(pd.read_csv(library[library_iterator], header = 0, delimiter="\t"))

    #search rowUnmapped
    for unmappedDF_iterator in range (rowUnmapped):
        found = False
        for library_list_iterator in range (len(librarydfs)):
            libRow, libCol = librarydfs[library_list_iterator].shape
            for library_df_iterator in range(libRow):
                if (unmappeddf.at[unmappedDF_iterator, "seq"]) in (librarydfs[library_list_iterator].at[library_df_iterator, "sequence"]).upper():
                    found = True
                    temp = pd.DataFrame([[librarydfs[library_list_iterator].at[library_df_iterator, "sgId"], unmappeddf.at[unmappedDF_iterator, "seq"], unmapped, unmappeddf.at[unmappedDF_iterator, "freq"], library[library_list_iterator]]], columns=output_columns)
                    outdf = outdf.append(temp)
        if (found == False and args.include_unfound):
            temp = pd.DataFrame([["N/A", unmappeddf.at[unmappedDF_iterator, "seq"], unmapped, unmappeddf.at[unmappedDF_iterator, "freq"], "N/A"]], columns=output_columns)
            outdf = outdf.append(temp)
        if args.percentage:
            sys.stdout.write('\r')
            sys.stdout.write("Percentage complete: " + str(unmappedDF_iterator/rowUnmapped * 100))

    #send output to csv
    outdf.to_csv(args.file_names[1], sep="\t", mode='a', index=False, header=header)

#call function
search_libraries(unmapped_file_name, library_file_names, header)
