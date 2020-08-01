#!/usr/bin/env python2.7
import sys

#load in nodes.dmp
ranks={}
for line in open(sys.argv[1]):
	cut=line.split('\t')
	ranks[cut[0]]=cut[4]

include=set(["species", "genus", "family", "order", "class", "phylum", "superkingdom"])

#prune blast output to remove mappings without a rank and remove taxid columnn
for  line in open(sys.argv[2]):
	cut=line.strip().split('\t')
	for i in range(len(cut)): 
		cut[i]=cut[i].strip()
	ids=cut[5]
	if len(ids.split(';'))<1: continue
	ct=0
	for id in ids.split(';'):
		if id not in ranks: continue
		if ranks[id] not in include: continue
		if ct>0: continue
		cut[5]=id
		ct+=1
		#print "\t".join(cut[:5] + cut[6:])
		#print "\t".join(cut[4:6])
		print "\t".join(cut)
