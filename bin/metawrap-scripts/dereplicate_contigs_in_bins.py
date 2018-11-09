#!/usr/bin/env python2.7
import sys, os

# Usage: ./script.py bins.stats binsFolder outFolder

# load in bin completion and contamination scores
print "Loading in bin completion and contamination scores..."
bin_scores={}
for line in open(sys.argv[1]):
	if "completeness" in line: continue
	cut=line.strip().split("\t")
	score = float(cut[1]) - 5*float(cut[2]) + 0.0000000001*float(cut[5])
	bin_scores[cut[0]]=score

# load in contigs in each bin
print "Loading in contigs in each bin..."
contig_mapping={}
for bin_file in os.listdir(sys.argv[2]):
	bin_name=".".join(bin_file.split("/")[-1].split(".")[:-1])
	for line in open(sys.argv[2]+"/"+bin_file):
		if line[0]!=">": continue
		contig=line[1:-1]
		if contig not in contig_mapping: contig_mapping[contig]=bin_name
		else: 
			if len(sys.argv)>4:
				if sys.argv[4]=="remove": contig_mapping[contig]=None
			elif bin_scores[bin_name]>bin_scores[contig_mapping[contig]]: 
				contig_mapping[contig]=bin_name


# go over the bin files again and make a new dereplicated version of each bin file	
print "Making a new dereplicated version of each bin file"
os.system("mkdir "+sys.argv[3])
for bin_file in os.listdir(sys.argv[2]):
        bin_name=".".join(bin_file.split("/")[-1].split(".")[:-1])
	out = open(sys.argv[3]+"/"+bin_file,'w')
	at_least_one=False
        for line in open(sys.argv[2]+"/"+bin_file):
		if line[0]==">":
			contig=line[1:-1]
			if contig_mapping[contig]==bin_name: 
				at_least_one=True
				store=True
			else: store=False
		if store==True: out.write(line)
	out.close()	
	if at_least_one==False:
		os.system("rm " + sys.argv[3] + "/" + bin_file)


