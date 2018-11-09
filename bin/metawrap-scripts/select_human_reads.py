#!/usr/bin/env python2.7
import sys

# This script takes in the reads that are proposed by bmtagger to be human, and filteres them out of the original fastq file.
# USAGE: ./select_human_reads.py bmreads.list reads_1.fastq > human_reads.fastq

# Load in the human reads:
human={}
for line in open(sys.argv[1]):
	human[line.strip()]=None

# only print the human reads
skip=True
for i, line in enumerate(open(sys.argv[2])):
	if i%4==0: 
		if line[1:].split("/")[0].split()[0] in human: skip=False
		else: skip=True 

	if skip==False: print line.rstrip()





