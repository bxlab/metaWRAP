#!/usr/bin/env python2.7
import sys, os


for filename in os.listdir(sys.argv[1]):
	bin_name=".".join(filename.split(".")[:-1])
	for line in open(sys.argv[1]+'/'+filename):
		if line[0]=='>':
			print line[1:-1]+'\t'+bin_name

