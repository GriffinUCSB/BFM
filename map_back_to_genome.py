#imports
import numpy as np
import pandas as pd
import argparse
import sys

#read gff3 as df
def process_gff3_to_df(filepath):
	file = open(r""+filepath,'r')
	txt = file.readlines()
	columns = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
	t = []
	for i in range(len(txt)):
		if ('#' not in txt[i]):
			t.append(txt[i].split('\t'))
	
	file.close()
	df = pd.DataFrame(data=t, columns = columns)
	return df

#read macs3 output into pandas df
def process_macs3_to_df(fp):
	file = open(r""+fp,'r')
	txt = file.readlines()
	columns = ['chr', 'start', 'end', 'length', 'abs_summit', 'pileup', '-log10(pvalue)', 'fold_enrichment', '-log10(qvalue)', 'name']
	t = []
	for i in range(len(txt)):
		if ('#' not in txt[i] and txt[i] != '\n'):
			t.append(txt[i].split('\t'))
	file.close()
	df = pd.DataFrame(data=t, columns = columns)
	df = df.drop([0], axis="index")
	df[["start", "end"]] = df[["start", "end"]].apply(pd.to_numeric)
	df = df.reset_index(drop=True)
	return df

#get gene name from info and keep transcripts only
def prepare_gff3(gff3):
	gff3 = gff3.loc[gff3['type'] == 'transcript']
	gff3 = gff3.reset_index(drop = True)
	gene_name = []
	for i in range(len(gff3)):
		name = gff3.at[i, 'attributes'].split('gene_name=')
		name = name[1].split(';')
		gene_name.append(name[0])
	gff3['gene'] = gene_name
	return gff3

#grab all the unique chromosomes in df
def create_chrom_list(df):
	e = pd.unique(df['chr']).tolist()
	return e

#subset the data
def subset(df, chr):
	#subset data
	subset = df.query('chr == @chr') #subset the data based on chromosome
	subset = subset.reset_index() #reset old index's
	#get feature start 
	subset = find_start_index(subset)
	return subset

#checks if feature is on top or bottom and gives an index to map to
def find_start_index(df):
	e = []
	for i in range(len(df)):
		if (df.at[i, 'strand'] == '+'):
			e.append(float(df.at[i,'start']))
		else:
			e.append(float(df.at[i,'end']))
	df['start_point'] = e
	return df


#calc midpoint
def midpoint(df):
	e = []
	for i in range(len(df)):
		e.append((df.at[i,'end'] + df.at[i,'start'])/2)
	df['midpoint'] = e
	return df

#find nearest function
def find_nearest(arr, value):
	arr = np.asarray(arr)
	d = np.abs(arr-value)
	idx = (d).argmin()
	return idx

#heavy lifter
def map(to_map,subsets):
	#empty arrays
	closest_feature = [] 
	distance = []
	#iterate through df
	for i in range(len(to_map)):
		#if chromosome in genome:
		if (to_map.at[i, 'chr'] in subsets[0]):

			#get index of correct chrmosome in df of chromosome subsets
			chrom_index = subsets[0].index(to_map.at[i, 'chr'])
			#get index of closest feature in chromosome df
			closest_feature_index = find_nearest(subsets[1][chrom_index][["start_point"]].to_numpy(), to_map.at[i, 'midpoint'])
			#take index and return gene symbol 
			closest_feature.append(subsets[1][chrom_index].at[closest_feature_index, 'gene'])
			#calculate distance between feature midpoint in chromosome and peak midpoint
			distance.append(np.abs(subsets[1][chrom_index].at[closest_feature_index, 'start_point'] -  to_map.at[i, 'midpoint'])/1000)
		#if chromosome not in geneome: skip peak
		else:
			closest_feature.append(" ")
			distance.append(" ")
	#add arrays to dataframe
	to_map['closest feature'] = closest_feature
	to_map['distance'] = distance
	#return dataframe
	return to_map


# argsparse functionality
parser = argparse.ArgumentParser(description='proccess user in')
parser.add_argument('--no_macs3', help="for files that are not macs3 output. File should be a .xlsx with at minimum columns labeled 'chr', 'start', and 'end'", action="store_true")
parser.add_argument('-g', '--gff3', nargs='+', default=[], help='input gff3 file (list of features by location)')
parser.add_argument('-p', '--peaks', nargs='+', default=[], help='input macs3 peaks output.')
args=parser.parse_args()
gff3_fp = args.gff3[0]
peaks_fp = args.peaks[0]
no_macs3 = args.no_macs3

# read in all features in the human genome
sys.stdout.write('Loading in gff3...')
sys.stdout.write('\n')
features = process_gff3_to_df(gff3_fp)

sys.stdout.write('Preparing gff3...')
sys.stdout.write('\n')
features = prepare_gff3(features)


# subset features based on chromosome
sys.stdout.write('Subsetting Data...')
sys.stdout.write('\n')
chr_list = create_chrom_list(features)
df_list = []
for i in chr_list:
	df_list.append(subset(features, i))
subsets = [chr_list, df_list]

#grab file to map peaks back to
sys.stdout.write('Loading in peaks...')
sys.stdout.write('\n')
if (no_macs3 == True):
	to_map = pd.read_excel(peaks_fp)
else:
	to_map = process_macs3_to_df(peaks_fp)

#create the annotation
sys.stdout.write('Mapping ...')
sys.stdout.write('\n')
to_map = midpoint(to_map)
to_map = map(to_map, subsets)
to_map = to_map.drop(columns=['midpoint'])

#output file
sys.stdout.write('Finished ...')
sys.stdout.write('\n')
output_temp = peaks_fp.split('.xls')
output_fp = output_temp[0] + '_annotated.txt'
to_map.to_csv(output_fp, sep="\t", mode='a', header=True, index=False)





