#!/usr/bin/env python2.7
import sys

# this script takes in a list of fastq read names and picks them out of a .fastq file.

#load in sequences of interest
seq={}
for line in open(sys.argv[1]):
	seq[line.strip()]=0

# only print out the fastq reads if they match one of the reads of interest
found=False
for i, line in enumerate(open(sys.argv[2])):
	line=line.strip()
	if i%4==0: 
		found=False
		read=line[1:].split("/")[0]
		if read in seq: found=True
		#for soi in seq: 
		#	if soi in line: found=True

	if found==True: 
		print line
		sys.stdout.flush()





