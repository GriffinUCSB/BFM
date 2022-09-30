import urllib.request
uniprot = "A6NKP2"

def uniprot_to_gene_name(uniprot):
    found = False
    with urllib.request.urlopen("http://www.uniprot.org/uniprot/" + uniprot + ".txt") as f:
    #print(data)
    #data.split("\n")
        i = 0
        while (found == False):
            line = f.readline()
            line = str(line)
            line.split('\n')
            if ("Name=" in line):
                #gene = line[4:]
                gene = line.split("Name=")[1]
                gene = gene.split(";")[0]
                found = True
            i+=1
    return(gene)
