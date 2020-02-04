#!/usr/bin/env python2.7
print "loading libs..."
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rc('font', family='arial')


def load_lib_sizes(filename):
	print "loading library sizes..."
	libs={}
	for line in open(filename):
		if line.startswith("#"): continue
		lib=line.strip().split("\t")[0]
		reads=int(line.strip().split("\t")[1])
		libs[lib] = reads
	return libs


def load_data(filename):
	print "loading abundance data..."
	df=pd.read_csv(filename, sep='\t', index_col=0)

	# remove all 0 rows
	df = df[(df.T != 0).any()]
	
	# standardize columns by total sum in each column
	df = df.div(df.sum(axis=0), axis=1)
	df=1000000*df
	return df


def set_colors_to_timeline(df):
	print "adding colored labels..."
	lut=[]
	for sample in df.columns.values:
		if "2013-04" in sample: lut.append('m')
        	elif "2014-09" in sample: lut.append('r')
        	elif "2015-06" in sample: lut.append('b')
		elif "2015-12" in sample: lut.append('c')
        	elif "2016-02" in sample: lut.append('y')
        	elif "2017-02" in sample: lut.append('g')
		else: lut.append('w')
	return lut


def draw_clustermap(df, lut):
	print "drawing clustermap..."
	sns.set(font_scale=1)
	df = df.fillna(0)
	if lut!=False:
		g = sns.clustermap(df, figsize=(14,8), col_colors=lut, col_cluster=True, yticklabels=True, cmap="magma")
	else:
		g = sns.clustermap(df, figsize=(14,8), col_cluster=True, yticklabels=True, cmap="magma")
	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
	plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)



# MAIN
df = load_data(sys.argv[1])

# log standardize:
df+=0.01; df=np.log(df)

# standardize rows by maximum value in each row
#df = df.div(df.max(axis=1), axis=0)


lut=False
draw_clustermap(df, lut)

if len(sys.argv)<3: out=".".join(sys.argv[1].split(".")[:-1])+".png"
else: out=sys.argv[2]


plt.savefig(out, bbox_inches='tight', dpi=300)


