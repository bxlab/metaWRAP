#!/usr/bin/env python2.7
import sys, os


# load bin contigs
contig_bins={}
for bin_file in os.listdir(sys.argv[1]):
	if bin_file.endswith(".fa") or bin_file.endswith(".fasta"): 
		bin_name=".".join(bin_file.split("/")[-1].split(".")[:-1])
		for line in open(sys.argv[1]+"/"+bin_file):
			if line[0]!=">": continue
			contig_bins[line[1:-1]]=bin_name

# store the read names and what bins they belong in in these dictionaries
# strict stores only perfectly aligning reads and permissive stores any aligned reads
mappings={}

# parse sam file from bwa mem -a alignment
for line in sys.stdin:
        if "NM:i:" not in line: continue
	read=line.split("\t")[0]
	contig=line.split("\t")[2]
	if contig not in contig_bins: continue
	if read not in mappings: mappings[read]=[{}, {}]
	if line.split("\t")[11]=="NM:i:0": mappings[read][0][contig_bins[contig]]=None
	mappings[read][1][contig_bins[contig]]=None

print mappings

