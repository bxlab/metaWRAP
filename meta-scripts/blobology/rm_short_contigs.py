#! /usr/bin/env python

import sys

f=open(sys.argv[1])
for line in f:
	if line.startswith(">"):
		cut=line.split("_")
		if int(cut[3])>999:
			print line.strip()
		else: quit
	else: print line.strip()
