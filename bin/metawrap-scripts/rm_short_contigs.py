#!/usr/bin/env python2.7
from __future__ import print_function
import sys


for line in open(sys.argv[2]):
	if not line.startswith(">"): print(line.strip())
	else:
		if int(line.split("_")[3])<int(sys.argv[1]): break
		else: print(line.strip())	
