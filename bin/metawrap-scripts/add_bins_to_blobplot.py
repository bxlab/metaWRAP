#!/usr/bin/env python2.7
# THis script takes in a blobplot text file and and number of bins as input, and annotates
# each line in the blob file with its bin name (if that contig exists in one of the bins)

import sys, os

# load the binned contigs:
contig_bins={}
for bin_file in os.listdir(sys.argv[2]):
	bin_name=bin_file.split("/")[-1]
	for line in open(sys.argv[2]+"/"+bin_file):
		if line.startswith(">"):
			contig_bins[line[1:-1]]=bin_name


# add bin annotations to blob file:
for line in open(sys.argv[1]):
	if line=="\n": continue
	if line.split("\t")[0]=="seqid":
		for i, field in enumerate(line.strip().split("\t")):
			if field=="taxlevel_phylum": phylum_column=i
		print line.strip() + "\tbin\tbinned_yes_no\tbinned_phylum"
	elif line.split("\t")[0] in contig_bins:
		phylum=line.split("\t")[phylum_column]
		print "\t".join([line.strip(), contig_bins[line.split("\t")[0]], "Binned", phylum])

	else:
		print "\t".join([line.strip(), "Unbinned", "Unbinned", "Unbinned"])

