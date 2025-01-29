#%% Initialize
import numpy as np
import pandas as pd
from datetime import datetime
now = datetime.now
import matplotlib.pyplot as plt
'''
This program uses a table of ATAC-Seq peaks (rows) and TFs (columns)
with 0s and 1s correesponding to whether the TF binds that ATAC-Seq peak
This program uses pre-clustered ATAC-Seq peaks (referred to as 'regions')
from Chromosome 1 and uses that to learn the clustering for regions from other chromosomes
'''
#%% readTable
df = pd.read_csv("AtacChip_OverlapTable.csv") # Input the table
clusters = []
table=[]
# add a row denoting the cluster identity of each 
for i in range(1,11): #10 Clusters
    f = open(f"10Clusters_Chr1/Cluster{i}.txt",'r')
    clusters.append([i.strip() for i in f.readlines()])
    f.close()
    tempdf = df[clusters[-1]].copy().T
    tempdf['cluster'] = i
    table.append(tempdf)
table = pd.concat(table)
table = table[['cluster']+list(table.columns)[:-1]]
clusters = [[i for row in clusters for i in row]]+clusters
del i, tempdf

#%% getProbabilities
# Get percentage of 1s for each TF in each cluster
lens = np.array([len(i) for i in clusters]) #Number of regions in each cluster
theta = np.zeros((table['cluster'].nunique(), len(table.columns)-1))
tfs = list(df.index)
for cluster in range(1,11):
    for tfi in range(len(tfs)):
        tf = tfs[tfi]
        vals = table[table['cluster']==cluster][tf]
        theta[cluster-1][tfi] = sum(vals)/len(vals)
del tf, tfi, vals, cluster

#%%  getOriginalClusterSize
size = []
for i in range(1,11):
    f = open(f"kMeansClustering/K562/10Clusters_Chr1/Cluster{i}.txt")
    size.append(len(f.readlines()))
    f.close()
del i, f

#%% clusterChromosome2
chrNo = 2
print(f"Clustering Chromosome {chrNo}")
if type(chrNo)==int: # Get regions of chromosome 2
    table2 = df[[i for i in df.columns if i[0]==chr(64+chrNo)]].T
else:
    table2 = df[[i for i in df.columns if i[0]==chrNo]].T
regions = list(table2.index)
table2['cluster'] = 0
l = len(table2)
for i in range(l):
    if i%1000 == l%1000: 
        print((l-i), now().time()) # Progress Tracker
    prob = [float(np.sum([np.log10(theta[j][t]) if table2.iloc[i,t]==1 else np.log10(1-theta[j][t])
                          for t in range(353)]))+np.log10(size[j]) for j in range(10)]
    # Get Maximum Likelihood Estimate for each region's clustering
    table2.iat[i, -1] = prob.index(max(prob))+1 #Assign the cluster to the region
table2.to_csv(f"kMeansClustering/K562/Learned Clusters/Learned_Chr{chrNo})withPrior.csv", header=True, index=True)
del i, prob
print(f"Clustered {chrNo}: {now().time()}")

#%% plotHeatmap
# Convert the table of 1s and 0s to an image for visualization
table2 = pd.read_csv(f"Learned Clusters/Learned_Chr{chrNo}.csv")
table2 = table2.set_index("Index")
atac = pd.read_csv("Raw Data/ATAC-Seq.csv")
atac = atac[atac['chr']==f'chr{chrNo}']
atac = atac.sort_values(by='start', ignore_index=True)
atacDict = atac[['name']].reset_index(drop=False).set_index('name').to_dict()['index']
val = table2[['cluster']].reset_index(drop=False)
val = val.sort_values(by=['cluster', 'Index'], key = np.vectorize(lambda i: int(i[1:-1]) if type(i)==str else i), ignore_index=True)
sortKey = val[['Index']].reset_index(drop=False).set_index('Index').to_dict()['index']
heatmap = table2.sort_index(key=np.vectorize(lambda i:sortKey[i]))
del heatmap['cluster']
img = heatmap.values.T
img = np.array([img[i//15] for i in range(len(img)*15)]) # Change aspect ratio
plt.imsave(f"Image Files/Learned from Chr1/Chromosome{chrNo}.png", img, cmap='binary') # black and white
print(f"Saved {chrNo}: {now().time()}")