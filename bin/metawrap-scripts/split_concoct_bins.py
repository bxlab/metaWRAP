#!/usr/bin/env python2.7
# Takes in the clustering_gt1000.csv file from CONCOCT binning, and splits the contigs into propper bins
# Usage:
# ./script clustering_gt1000.csv assembly_file.fa out_folder

import sys, os

print "Loading in the bins that the contigs belong to..."
bins={}
for line in open(sys.argv[1]):
	if line.startswith("contig_id"):
		continue
	bins[line.strip().split(",")[0].split(".")[0]] = line.strip().split(",")[1]


print "Going through the entire assembly and splitting contigs into their respective bin file..."
current_bin=""
for line in open(sys.argv[2]):
	if line.startswith(">"):
		if current_bin!="": 
			f.close()
		contig = line[1:-1].split(".")[0].split()[0]
		line=line.rsplit()[0]+"\n"
		if contig in bins: 
			current_bin="bin."+bins[contig]+".fa"
		else: 
 			current_bin="unbinned.fa"
		f = open(sys.argv[3]+"/"+current_bin,'a')
	f.write(line)
print "Done!"
	
	

		
