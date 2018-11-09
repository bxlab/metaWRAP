#!/usr/bin/env python2.7
import sys

shorten=False
for line in open(sys.argv[1]):
	if line[0]!=">":
		print line.rstrip()
	else:
		if shorten==True:
			print "_".join(line.rstrip().split("_")[:4])
		elif len(line)>20 and len(line.split("_"))>5:
			print "_".join(line.rstrip().split("_")[:4])
			shorten=True
		else:
			print line.rstrip()


