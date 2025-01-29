import pandas as pd
import numpy as np

'''
This program loads a csv file containing Protein IDs as one of the columns,
reads the corresponding Protein Sequence from another file, calculates its
Shannon Entropy and adds it to another column in the csv
'''
def get_entropy(seq):
	d = {chr(i):0 for i in range(65,65+26)}
	for i in seq:
		d[i]+=1
	p = np.array([i for i in d.values() if i > 0])
	p = p/sum(p)
	logp = np.log2(p)
	return -sum(p*logp)


d = {} #Dictionary of SeqID to Seq
f = open("Raw Data/Over60_ChromosomeAll.fasta").readlines() #File with Seqs
name = ""; seq = ""
for i in f:
	if i[0]==">":
		d[name]=seq
		name = i.strip()[1:]
		seq = ""
	else:
		seq+=i.strip()
d[name]=seq

df = pd.read_csv('TopHitsOnly_noCutoff.csv',sep="$") #csv File
df['Seq'] = df['Query'].apply(lambda i: d[i])
df['entropy'] = df['Seq'].apply(get_entropy)
df.to_csv('TopHitsOnly_noCutoff.csv',sep="$", header = True, index = False)