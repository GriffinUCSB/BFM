#find unmapped reads V2 6/14/2021
#Griffin Kramer

#!important note! the built in delimiter for unmapped is a space based on the samtools output of the line:
# samtools view -F 4 file.bam | awk '{print $10}' | sort | uniq -c | sort -k1nr | head -n 20 > unmapped.txt
#the delimiter for librarys is a tab ("\t")
#both of these can easily be changed manually below

#imports
import csv, pandas as pd
import os

#variables for configuration
unmapped_file_name = "fixedfinalUnmapped.txt"
library_file_names = ["urLib1.txt", "urLib2.txt"]
header = True #If this program is being used to create a new .csv, set header to True in config

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

    #search
    for unmappedDF_iterator in range (rowUnmapped):
        found = False
        for library_list_iterator in range (len(librarydfs)):
            libRow, libCol = librarydfs[library_list_iterator].shape
            for library_df_iterator in range(libRow):
                if (unmappeddf.at[unmappedDF_iterator, "seq"]) in (librarydfs[library_list_iterator].at[library_df_iterator, "sequence"]).upper():
                    found = True
                    temp = pd.DataFrame([[librarydfs[library_list_iterator].at[library_df_iterator, "sgId"], unmappeddf.at[unmappedDF_iterator, "seq"], unmapped, unmappeddf.at[unmappedDF_iterator, "freq"], library[library_list_iterator]]], columns=output_columns)
                    outdf = outdf.append(temp)
        #if (found == False):
        #    temp = pd.DataFrame([["N/A", unmappeddf.at[unmappedDF_iterator, "seq"], unmapped, unmappeddf.at[unmappedDF_iterator, "freq"], "N/A"]], columns=output_columns)
        #    outdf = outdf.append(temp)

    #send output to csv
    outdf.to_csv('foundUnmapped3.txt', sep="\t", mode='a', index=False, header=header)

#call function
search_libraries(unmapped_file_name, library_file_names, header)
