import pandas as pd
import sys

proteomics = "20220830xDonorMSstats.xlsx"
mapped_uniprot = "uniprot_remapped.xlsx"
data = ["ds_ZAPvGFP_MAGeCK-iNC_Output/pref_ds_ZAPvGFP_all_genes.csv","ss_PREvGFP_MAGeCK-iNC_Output/pref_ss_PREvGFP_all_genes.csv"]
header = True

data1df = pd.read_csv(data[0])
data2df = pd.read_csv(data[1])

prot_df = pd.read_excel(proteomics, index_col = 0)
prot_df = prot_df.reset_index(drop=True)
mapped_df = pd.read_excel(mapped_uniprot, usecols = "A,C")
mapped_df.columns = ["Protein", "Gene"]

temp = []
#for i in range():
for i in range (len(mapped_df)):
    first = True
    for j in range (len(data1df["gene"])):
        if (data1df.at[j, "gene"] in mapped_df.at[i, "Gene"]) and (first == True):
            #mapped_df.at[i, "Gene"] = data1df.at[j, "gene"]
            temp.append(data1df.at[j, "product"])
            first = False
            sys.stdout.write('\r')
            sys.stdout.write("1: " +  str(i) + "/" + str(len(mapped_df)))
    if first == True:
        temp.append(" ")


mapped_df["DS_Product"] = temp

temp = []
for i in range (len(mapped_df)):
    first = True
    for j in range (len(data2df["gene"])):
        if (data2df.at[j, "gene"] in mapped_df.at[i, "Gene"]) and (first == True):
            #mapped_df.at[i, "Gene"] = data2df.at[j, "gene"]
            temp.append(data2df.at[j, "product"])
            first = False
            sys.stdout.write('\r')
            sys.stdout.write("2: " +  str(i) + "/" + str(len(mapped_df)))
    if first == True:
        temp.append(" ")


mapped_df["SS_Product"] = temp


prot_df = prot_df.merge(mapped_df, on="Protein")

prot_df.to_excel("20220830xDonorMSstats_annotated.xlsx")
