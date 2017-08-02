#!/usr/bin/env python
# THis script takes in a blobplot text file and and number of bins as input, and annotates
# each line in the blob file with its bin name (if that contig exists in one of the bins)

import sys

# load the binned contigs:
contig_bins={}
for bin in sys.argv[2:]:
	bin_name=bin.split("/")[-1]
	for line in open(bin):
		if line.startswith(">"):
			contig_bins[line[1:-1]]=bin_name


# add bin annotations to blob file:
for line in open(sys.argv[1]):
	if line.split("\t")[0]=="seqid":
		print line.strip() + "\t" + "bins"
	elif line.split("\t")[0] in contig_bins:
		print line.strip() + "\t" + contig_bins[line.split("\t")[0]]
	else:
		cut=line.strip().split("\t")
		for i in range(4, len(cut)):
			cut[i] = "Not annotated"
		cut.append("Not annotated")
		print "\t".join(cut)

