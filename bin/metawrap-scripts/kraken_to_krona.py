#!/usr/bin/env python2.7
import sys
# This script takes in a translated kraken file of either contigs (from SPAdes) or reads, and parses it into a format for ktImportText to produce a kronachart.
data={}
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	if len(cut)<2:
		continue
	name=cut[0]
	tax="\t".join(cut[1].split(";"))
	
	#this is a contig! calculate weight with cov*length
	if "length" in name and len(name.split("_"))>5:
		if tax in data: data[tax]+=float(name.split("_")[3]) * float(name.split("_")[5])
                else: data[tax]=float(name.split("_")[3]) * float(name.split("_")[5])

	
	#this is a read! weight is 1
	else:
		if tax in data: data[tax]+=1
		else: data[tax]=1




for tax in data:
	print str(data[tax]) + "\t" + tax

