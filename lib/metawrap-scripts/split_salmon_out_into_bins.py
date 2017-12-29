#!/usr/bin/env python
import sys, os

# load bins and their contigs
bins={}
for bin_name in os.listdir(sys.argv[2]):
	for line in open(sys.argv[2]+"/"+bin_name):
		if line[0]!=">": continue
		bins[line[1:-1]]=bin_name

# load in the salmon abundances
bin_abundances={}
for salmon_file in os.listdir(sys.argv[1]):
	sample=salmon_file.split(".")[0]
	for line in open(sys.argv[1]+'/'+salmon_file):
		if "transcript" in line: continue
		if line.split("\t")[0] not in bins: continue
		bin=bins[line.split("\t")[0]]
		abun=float(line.strip().split("\t")[1])
	
		if bin not in bin_abundances:
			bin_abundances[bin]={}
			bin_abundances[bin]["total_len"]=0
			bin_abundances[bin]["total_cov"]=0
			bin_abundances[bin]["samples"]={}
	
		bin_abundances[bin]["total_len"]+=int(line.split("\t")[0].split("_")[3])
		bin_abundances[bin]["total_cov"]+=abun * int(line.split("\t")[0].split("_")[3])
	
	
	for bin in bin_abundances:
		bin_abundances[bin]["samples"][sample]=bin_abundances[bin]["total_cov"] / bin_abundances[bin]["total_len"]	
                bin_abundances[bin]["total_len"]=0
                bin_abundances[bin]["total_cov"]=0


first=True
for bin in bin_abundances: 
	if first==True:
		sys.stdout.write('Genomic bins')
		for sample in bin_abundances[bin]["samples"]:
			sys.stdout.write('\t'+sample)
		print ""
		first=False
	sys.stdout.write('.'.join(bin.split('.')[:-1]))
	for sample in bin_abundances[bin]["samples"]:
		sys.stdout.write("\t" + str(bin_abundances[bin]["samples"][sample]))
	print ""




