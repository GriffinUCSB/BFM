import pandas as pd
import uniprotTOgenesymbol as utgs
import sys
proteomics = "20220830xDonorMSstats.xlsx"
prot_df = pd.read_excel(proteomics, index_col = 0)
outdf = prot_df["Protein"]
outdf.to_csv("uniprot_list.txt", sep="\t", mode='a', index=False,header = True)
