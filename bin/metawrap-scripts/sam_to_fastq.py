#!/usr/bin/env python2.7
import sys
for line in open(sys.argv[1]):
	cut=line.split('\t')
	if len(cut)<11: continue
	print ">"+cut[0]
	print cut[9]
	print "+"
	print cut[10]
