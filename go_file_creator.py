import pandas as pd
import sys 

#Get id per genesymbol
Go_id = open(r"goa_human.gaf", 'r')
Go_id_txt = Go_id.readlines()
id_columns  = ["Database", "Object ID", "Object Symbol", "Qualifier", "GO ID", "Refrence", "Evidence Code", "With/From", "Aspect"]
t = []
data = []
c = 0
for i in range (len(Go_id_txt)):
	z = []
	if ('!' not in Go_id_txt[i]):
		t.append(Go_id_txt[i].split('\t'))
		for j in range(9):
			z.append(t[c][j])
		data.append(z)
		c += 1
Go_id.close()
Go_id_df = pd.DataFrame(data = data, columns = id_columns)


#get name per id
Go_name = open(r"go-basic.obo", 'r')
Go_name_txt = Go_name.readlines()
columns = ['id', 'name', 'namespace']
t = []
for i in range (len(Go_name_txt)):
	if ('[Term]' in Go_name_txt[i]):
		temp_id = Go_name_txt[i+1].split('\n')[0].split('id: ')[1]
		temp_name =  Go_name_txt[i+2].split('\n')[0].split('name: ')[1]
		temp_namespace = Go_name_txt[i+3].split('\n')[0].split('namespace: ')[1]
		t.append([temp_id, temp_name, temp_namespace])
Go_name_df = pd.DataFrame(data = t, columns = columns)

def subset(df, this_id):
	subset = df.query('id == @this_id')
	subset = subset.reset_index() #reset old index's
	return subset.at[0, "name"], subset.at[0, 'namespace']

#Create new file
Go_name = []
GO_namespace = []
for i in range (len(Go_id_df)):
	sys.stdout.write('\r')
	sys.stdout.write('mapping ' + str(i/len(Go_id_df)))
	this_id = Go_id_df.at[i, "GO ID"]
	t_name, t_namespace = subset(Go_name_df, this_id)
	Go_name.append(t_name)
	GO_namespace.append(t_namespace)
Go_id_df['GO Name'] = Go_name
Go_id_df['GO_Namespace'] = GO_namespace
Go_id_df.to_csv('go_super.txt', sep="\t", mode='a', header=True, index=False)
#Go_supreme = 
