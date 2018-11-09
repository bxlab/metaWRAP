#!/usr/bin/env python2.7
import sys

# This script takes in the reads that are probosed by bmtagger to be human, and filteres them out of the original fastq file.
# USAGE: ./skip_human_reads.py bmreads.list reads_1.fastq > reads_1.clean.fastq

# Load in the human reads:
human={}
for line in open(sys.argv[1]):
        human[line.strip()]=None

# Print out the fastq file line by line unless the read is human:
skip=False
for i, line in enumerate(open(sys.argv[2])):
        if i%4==0:
                if line[1:].split("/")[0].split()[0] in human: skip=True
		else: skip=False

        if skip==False: print line.rstrip()

