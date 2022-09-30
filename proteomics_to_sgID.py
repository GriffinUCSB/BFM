#to use script, delimiters and index of start may need to be changed -> see below
import pandas as pd
import uniprotTOgenesymbol as utgs
import sys
#argsparse stuff incoming
raw = "PooledScreenGuides.csv"
lib = "20150624_CRISPRiv2_final.Top5.txt"
header = True

def check_gs(gene, uniprot):
    if pd.isna(gene):
        gene = utgs.uniprot_to_gene_name(uniprot)
        return gene
        #print("uniprot check")
    else:
        return gene
    #if ";" in gene:
        #print('hello')

def process_gs(raw, lib, header):
#load files
    # "," delimiter for csv, "\t for .txt files"
    tbl = pd.read_csv(raw, header=None, delimiter=",")
    tbl.columns = ['uniprot', 'symbol', 'info']
    rows, cols = tbl.shape
    #index of start. look at protoeomics data file output, set to start of data
    tbl = tbl.iloc[7:]
    #print(tbl)

    lib = pd.read_csv(lib, header=None, delimiter='\t')
    lib.columns = ['sgID', 'flag2', 'flag3', 'space', 'gene', 'flag4', 'sequence', 'flag7', 'complement', 'flag9']
    #print(lib)
    #prepare output
    output_cols = ['sgID', 'sequence', 'gene']
    outdf = pd.DataFrame(columns = output_cols)

    #print (tbl)
    #print (tbl.at[25, 'symbol'])
    for i in range(rows-7):
        gene = tbl.at[i+7, 'symbol']
        gene = check_gs(gene, tbl.at[i+7,'uniprot'])
        if ";" in gene:
            genelist = gene.split(";")
            gene1 = genelist[0]
            gene2 = genelist[1]
            out1 = lib.query("gene == @gene1")
            out2 = lib.query("gene == @gene2")
            out1 = out1.reset_index()
            out2 = out2.reset_index()
            out = pd.concat([out1,out2], ignore_index = True)
        #print(str(i) + "/" + str(rows))
        #print (gene)
        else:
            out = lib.query("gene == @gene")
            out = out.reset_index()
        #print(out)
        for j in range(len(out)):
            temp = pd.DataFrame([[out.at[j, 'sgID'], out.at[j,'sequence'], out.at[j,'gene']]], columns=outdf.columns)
            #print(temp)
            outdf = pd.concat([outdf,temp], ignore_index = True)

        sys.stdout.write('\r')
        sys.stdout.write("Percentage complete: " + str(i/rows*100)[:4])

    outdf.to_csv("sgID_list.txt", sep="\t", mode='a', index=False, header=header)




process_gs(raw, lib, header)
