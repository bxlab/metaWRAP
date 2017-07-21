#!/usr/bin/env python
import sys
import textwrap

dic={}
tmp_contig=""
for line in open(sys.argv[1]):
	if line[0]==">": 
		if tmp_contig!="":
			dic[name]=tmp_contig
			tmp_contig=""
		cut=line.strip()[1:].split(" ")
		name=">"+cut[0]+"_length_"+cut[3].split("=")[1]+"_cov_"+cut[2].split("=")[1]
	else: tmp_contig+=line.strip()	
dic[name]=tmp_contig


for k in sorted(dic, key=lambda k: len(dic[k]), reverse=True):
        print k
        print textwrap.fill(dic[k], 100, break_on_hyphens = False)
