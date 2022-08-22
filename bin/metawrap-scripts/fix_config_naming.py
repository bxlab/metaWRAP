#!/usr/bin/env python2.7
from __future__ import print_function
import sys

for line in open(sys.argv[1]):
	if line[0]==">":
		for c in line:
			if c=="=": c="_"
			sys.stdout.write(c)
	else:
		print(line.rstrip())
