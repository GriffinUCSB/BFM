#imports
import csv, pandas as pd
import os

#variables for configuration
count_file_name = "ds1A.sort.txt"
library_file_names = ["urLib2.txt", "urLib1.txt"]
header = True #If this program is being used to create a new .csv, set header to True in config

def create_custom(count, library, header):

    #prepare output dataframe
    output_columns = ["sgId", "sublibrary", "gene", "transcripts", "sequence"]
    outdf = pd.DataFrame(columns = output_columns)

    failed_columns = ["sgId"]
    faileddf = pd.DataFrame(columns = failed_columns)

    #load count table .txt file into dataframe
    counttbl = pd.read_csv(count, header = None, delimiter="\t")
    counttbl.columns = ['feature', 'freq']
    rows, cols = counttbl.shape

    #load libraries into dataframes
    librarydfs = []
    for library_iterator in range (len(library)):
        librarydfs.append(pd.read_csv(library[library_iterator], header = 0, delimiter="\t"))

    #do the search
    for count_i in range (rows):
        found = False
        print(str(count_i) + "/" + str(rows))

        for library_i in range (len(librarydfs)):

            libRow, libCol = librarydfs[library_i].shape
            library_df_i = 0

            while (library_df_i < libRow and found == False):

                if (counttbl.at[count_i, "feature"]) in (librarydfs[library_i].at[library_df_i, "sgId"]):
                    found = True
                    temp = pd.DataFrame([[counttbl.at[count_i, "feature"], "custom", librarydfs[library_i].at[library_df_i, "gene"], librarydfs[library_i].at[library_df_i, "transcripts"], librarydfs[library_i].at[library_df_i, "sequence"]]], columns=output_columns)
                    outdf = outdf.append(temp)

                library_df_i+=1

        if (found == False):
            temp2 = pd.DataFrame([[counttbl.at[count_i, "feature"]]], columns = failed_columns)
            faileddf = faileddf.append(temp2)

    outdf.to_csv('custom_lib_table.txt', sep="\t", mode='a', index=False, header=header)
    faileddf.to_csv('unfound_features_table.txt', sep="\t", mode='a', index=False, header=header)

#call function
create_custom(count_file_name, library_file_names, header)
