#the purpose of this script is to create a MAGeCK INC compatible count table through merging
#imports
import csv, pandas as pd
import os

#input count file names here
count_file_names = ["B1-T0.sort.txt", "B1-Un-T14.sort.txt", "B1-Z-T14.sort.txt", "B2-T0.sort.txt", "B2-Un-T14.sort.txt", "B2-Z-T14.sort.txt"]
header = True #If this program is being used to create a new .csv, set header to True in config

def get_gene_name(sgID):
    if ("neg" in sgID):
        out = "non-targeting"
    elif ("_+_" in sgID):
        out = sgID.split('_+_')[0]
    else:
        out = sgID.split('_-_')[0]
    return out

#Step 1: merge some mf count tables
def create_wide_count_table(names, header):

    #create header for output file
    output_columns = ["sgRNA", "gene"]
    #output_columns = ["sgRNA", gene]
    for i in range (len(names)):
        output_columns.append(names[i])
    outdf = pd.DataFrame(columns = output_columns)

    # load in first table to get length
    tbl = pd.read_csv(names[0], header = None, delimiter="\t")
    tbl.columns = ['sgRNA', names[0]]
    #tbl.at[names[0]] = tbl.at['genes']
    rows, cols = tbl.shape


    for k in range(len(names)-1):
        temp_df = pd.read_csv(names[k+1], header = 0, delimiter="\t")
        temp_df.columns = ['sgRNA', names[k+1]]
        tbl = tbl.merge(temp_df, on='sgRNA')
        print(str(k)+"/"+str(len(names)))

    temp_array = []
    for i in range(rows-1):
        temp_array.append(get_gene_name(tbl.at[i, "sgRNA"]))
    #print(temp_array)
    tbl.insert(1, "gene", pd.Series(temp_array))
    #for i in range(rows-1):
    #    tbl.at[i, "gene"] = get_gene_name(tbl.at[i, "sgRNA"])

    outdf = outdf.append(tbl)
    outdf.to_csv('PEX_Screen_wct.txt', sep="\t", mode='a', index=False, header=header)

create_wide_count_table(count_file_names, header)
