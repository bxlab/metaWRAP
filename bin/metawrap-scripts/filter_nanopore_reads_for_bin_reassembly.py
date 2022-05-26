#!/usr/bin/env python2.7
import sys, os

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n', '*':'*'} 
def rev_comp(seq):
	rev_comp=""
	for n in seq:
		rev_comp+=complement[n]
	return rev_comp[::-1]

# load bin contigs
print "loading contig to bin mappings..."
contig_bins={}
for bin_file in os.listdir(sys.argv[1]):
	if bin_file.endswith(".fa") or bin_file.endswith(".fasta"): 
		bin_name=".".join(bin_file.split("/")[-1].split(".")[:-1])
		for line in open(sys.argv[1]+"/"+bin_file):
			if line[0]!=">": continue
			contig_bins[line[1:-1]]=bin_name

# store the read names and what bins they belong in in these dictionaries
# strict stores only perfectly aligning reads and permissive stores any aligned reads

print "Parsing sam file and writing reads to appropriate files depending what bin they alligned to..."
files={}
opened_bins={}
for line in sys.stdin:
	if line[0]=="@": continue
	cut = line.strip().split("\t")
	binary_flag = bin(int(cut[1]))

	# skip non aligned reads
	if cut[2]=="*": continue

	# make sure the R and F reads aligned to the same bin
	if cut[2] not in contig_bins: continue

	# make sure the reads aligned again
	if "NM:i:" not in line: continue

	bin_name = contig_bins[cut[2]]
	# open the revelant output files
	if bin_name not in opened_bins:
		opened_bins[bin_name]=None
		files[sys.argv[2]+"/"+bin_name+".nanopore.fastq"]=open(sys.argv[2]+"/"+bin_name+".nanopore.fastq", "w")

	# determine alignment type from bitwise FLAG
	binary_flag = bin(int(cut[1]))

	# if the reads are reversed, fix them
	try:
		if binary_flag[-5]=='1':
			cut[9] = rev_comp(cut[9])
			cut[10] = cut[10][::-1]
	except IndexError:
		pass


	# strict assembly
	files[sys.argv[2]+"/"+bin_name+".nanopore.fastq"].write('@' + cut[0] + "/1" + "\n" + cut[9] + "\n+\n" + cut[10] + "\n")


print "closing files"
for f in files:
	files[f].close()


print "Finished splitting reads!"
