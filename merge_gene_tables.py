import csv, pandas as pd
import os

file = "PexScreen_genetable_EXL.xls"
file2 = "ssds_combined_genetable_EXL.xls"
header = True

def merge(og_file, file_to_merge, header):

    table = pd.read_excel(og_file, header = [0,1,2], index_col = 0)
    table2 = pd.read_excel(file_to_merge, header = [0,1,2], index_col = 0)
    rows,cols = table2.shape

    gene_names_arr=table.index.values
    gene_names_arr2 = table2.index.values
    ctr = 0

    #print(table.iloc[3])

    tempdf=pd.DataFrame(columns = table2.columns)
    for i in range (len(gene_names_arr)):
        if (ctr >= len(gene_names_arr2)):
            na_arr = []
            for j in range(cols):
                na_arr.append(0)
            tempdf2 = pd.DataFrame([na_arr], columns = table2.columns, index=[gene_names_arr[i]])
            tempdf = tempdf.append(tempdf2)
        elif (gene_names_arr[i] in gene_names_arr2[ctr]):
            tempdf = tempdf.append(table2.iloc[ctr])
            ctr = ctr+1
        else:
            na_arr = []
            for j in range(cols):
                na_arr.append(0)
            tempdf2 = pd.DataFrame([na_arr], columns = table2.columns, index=[gene_names_arr[i]])
            tempdf = tempdf.append(tempdf2)
        print(str(i) + "/" + str(len(gene_names_arr)))
    #print(tempdf)


    table = table.join(tempdf)
    table.to_csv('supermerge_genetable.txt', sep="\t", mode='a', index=True, header=header)


merge(file, file2, header)
