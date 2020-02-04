#!/usr/bin/env python2.7
import sys, os
import numpy as np

# load bins and their contigs
bins={}
for bin_name in os.listdir(sys.argv[2]):
	for line in open(sys.argv[2]+"/"+bin_name):
		if line[0]!=">": continue
		bins[line[1:-1]]=bin_name

# load contig lengths:
contig_lengths={}
length=0
for line in open(sys.argv[3]):
	if line[0]==">":
		if length!=0: 
			contig_lengths[contig]=length
		length=0
		contig=line[1:-1].split()[0]
	else:
		length+=len(line.strip())
contig_lengths[contig]=length


# load in the salmon abundances
bin_abundances={}
for salmon_file in os.listdir(sys.argv[1]):
	sample=".".join(salmon_file.split(".")[:-2])
	ct=0
	for line in open(sys.argv[1]+'/'+salmon_file):
		if "transcript" in line: continue
		if line.split("\t")[0] not in bins: continue
		ct+=1
		bin=bins[line.split("\t")[0]]
		abun=float(line.strip().split("\t")[1])
		if bin not in bin_abundances:
			bin_abundances[bin]={}
			bin_abundances[bin]["total_len"]=0
			bin_abundances[bin]["total_cov"]=0
			bin_abundances[bin]["cov_list"]=[]
			bin_abundances[bin]["samples"]={}


		length = contig_lengths[line.split("\t")[0]]
		weight = length/1000

		for i in range(weight):
			bin_abundances[bin]["cov_list"].append(abun)
		bin_abundances[bin]["total_len"]+=contig_lengths[line.split("\t")[0]]
		bin_abundances[bin]["total_cov"]+=abun * contig_lengths[line.split("\t")[0]]
	
	if ct==0:
		sys.stderr.write("\nNone of the contigs/scaffolds in the -a metagenomic assembly file were present in the bin files. Please make sure that the bins and total assembly have the exact same bins. One cause for this could be that you reassembled the bins, disrupring the contig naming. If you do not have the original total metagenomic assembly file, then you could not provide the -a option at all (but this is not ideal for abundance estimation).\n")
		sys.exit(1)
		quit()

	for bin in bin_abundances:
		#bin_abundances[bin]["samples"][sample]=bin_abundances[bin]["total_cov"] / bin_abundances[bin]["total_len"]	
		bin_abundances[bin]["samples"][sample] = np.median(bin_abundances[bin]["cov_list"])
		#bin_abundances[bin]["samples"][sample] = sum(bin_abundances[bin]["cov_list"])
                bin_abundances[bin]["total_len"]=0
                bin_abundances[bin]["total_cov"]=0
		bin_abundances[bin]["cov_list"]=[]


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




