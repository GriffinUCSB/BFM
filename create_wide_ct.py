#imports
import csv, pandas as pd
import os

#input count file names here
count_file_names = ["ds1A.sort.txt","ds1B.sort.txt","ds1C.sort.txt","ds1D.sort.txt","ds1E.sort.txt","ds2A.sort.txt","ds2B.sort.txt","ds2C.sort.txt","ds2D.sort.txt",
"ds2E.sort.txt", "ss1llib.sort.txt", "ss1lbfp.sort.txt","ss1lgfp.sort.txt","ss1lpre.sort.txt",
"ss1lnon.sort.txt","ss2llib.sort.txt", "ss2lbfp.sort.txt","ss2lgfp.sort.txt","ss2lpre.sort.txt","ss2lnon.sort.txt",
"sskllib.sort.txt", "ssklbfp.sort.txt","ssklgfp.sort.txt","ssklpre.sort.txt","ssklnon.sort.txt"]
header = True #If this program is being used to create a new .csv, set header to True in config

def get_gene_name(sgID):
    if ("neg" in sgID):
        out = "non-targeting"
    elif ("_+_" in sgID):
        out = sgID.split('_+_')[0]
    else:
        out = sgID.split('_-_')[0]
    return out

def create_wide_count_table(names, header):

    #create header for output file
    output_columns = ["sgRNA", "gene"]
    for i in range (len(names)):
        output_columns.append(names[i])
    outdf = pd.DataFrame(columns = output_columns)

    # load in first table to get length
    first_tbl = pd.read_csv(count_file_names[0], header = None, delimiter="\t")
    rows, cols = first_tbl.shape

    #do the job
    for j in range(rows-1):

        first = True
        temp_df2 = pd.DataFrame(columns = output_columns)

        for k in range(len(names)):

            temp_df = pd.read_csv(names[k], header = 0, delimiter="\t")
            temp_df.columns = ['sgRNA', 'freq']

            if (first == True):
                temp_df2.at[j, "sgRNA"] = temp_df.at[j, "sgRNA"]
                temp_df2.at[j, "gene"] = get_gene_name(temp_df2.at[j, "sgRNA"])
                first = False

            temp_df2.at[j, names[k]] = temp_df.at[j, "freq"]

        #print(temp_df2)
        outdf = outdf.append(temp_df2)
        print(str(j)+"/"+str(rows))

    outdf.to_csv('wide_count_table_FIXED.txt', sep="\t", mode='a', index=False, header=header)

create_wide_count_table(count_file_names, header)
