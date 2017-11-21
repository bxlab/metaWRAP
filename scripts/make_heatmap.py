#!/usr/bin/env python


import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# load in library sizes:
libs={}
for line in open(sys.argv[1]):
	if line.startswith("#"): continue
	lib=line.strip().split("\t")[0]
	reads=int(line.strip().split("\t")[1])
	libs[lib] = reads
	

# load in abundannce table as pandas dataframe
z=pd.read_csv(sys.argv[2], sep='\t', index_col=0)
z.columns.name

# remove all 0 rows
z = z[(z.T != 0).any()]

# standardize abundance values by total number of reads in each sample
for column in z: z[column]=1000000*z[column]/libs[column]


# standardize columns by total sum in each column
z = z.div(z.sum(axis=0), axis=1)
z=1000000*z

# standardize rows by maximum value in each row
#z = z.div(z.max(axis=1), axis=0)

# add colored x-labels
lut=[]
for sample in z.columns.values:
	#skips SG1-N samples:
	if "-N" in sample:
		z = z.drop(sample, 1)
		continue

	if "2013-04" in sample: lut.append('m')
        elif "2014-09" in sample: lut.append('r')
        elif "2015-06" in sample: lut.append('b')
	elif "2015-12" in sample: lut.append('c')
        elif "2016-02" in sample: lut.append('y')
        elif "2017-02" in sample: lut.append('g')
	else: lut.append('w')


# log standardize:
z=z+0.01
z = np.log(z)

# make heat map
sns.set(font_scale=1)
g = sns.clustermap(z, figsize=(14,8), col_colors=lut, col_cluster=True)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)


plt.savefig(sys.argv[3], bbox_inches='tight')


