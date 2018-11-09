#!/usr/bin/env python2.7
import sys
import textwrap

dic={}
tmp_contig=""
for line in open(sys.argv[1]):
	if line[0]==">": 
		if tmp_contig!="":
			dic[name]=tmp_contig
			tmp_contig=""
		name=line.strip()
	else: tmp_contig+=line.strip()	
dic[name]=tmp_contig


for k in sorted(dic, key=lambda k: len(dic[k]), reverse=True):
        print k
	print textwrap.fill(dic[k], 100, break_on_hyphens = False)
