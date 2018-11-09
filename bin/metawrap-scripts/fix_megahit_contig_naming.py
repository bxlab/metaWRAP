#!/usr/bin/env python2.7
import sys
import textwrap

dic={}
tmp_contig=""
min_len=int(sys.argv[2])
good=True
for line in open(sys.argv[1]):
	if line[0]==">": 
		if tmp_contig!="":
			if good==True: 
				dic[name]=tmp_contig
			tmp_contig=""
		cut=line.strip()[1:].split(" ")
		
		if int(cut[3].split("=")[1])<min_len:
			good=False
		else: 
			good=True
		
		name=">"+cut[0]+"_length_"+cut[3].split("=")[1]+"_cov_"+cut[2].split("=")[1]
	else: tmp_contig+=line.strip()	
dic[name]=tmp_contig


for k in sorted(dic, key=lambda k: len(dic[k]), reverse=True):
        print k
        print textwrap.fill(dic[k], 100, break_on_hyphens = False)
