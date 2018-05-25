#!/usr/bin/env python
#usage: 
# bwa mem -a assembly.fa reads_1.fastq reads_2.fastq | ./filter_reads_for_bin_reassembly.py original_bin_folder reads_1.fastq reads_2.fastq output_dir
import sys, os
cache = 100

# load bin contigs
contig_bins={}
for bin_file in os.listdir(sys.argv[1]):
	if bin_file.endswith(".fa") or bin_file.endswith(".fasta"): 
		bin_name=".".join(bin_file.split("/")[-1].split(".")[:-1])
		for line in open(sys.argv[1]+"/"+bin_file):
			if line[0]!=">": continue
			contig_bins[line[1:-1]]=bin_name

# store the read names and what bins they belong in in these dictionaries
# strict stores only perfectly aligning reads and permissive stores any aligned reads
mapping={}


print "\nParsing sam file from bwa mem -a alignment. Warning: the IDs of all the reads mapped to the bins will be stored in memory (might take up quite a bit of RAM)"
for line in sys.stdin:
        if "NM:i:" not in line: continue
	read=line.split("\t")[0]
	contig=line.split("\t")[2]
	if contig not in contig_bins: continue
	if read not in mapping: mapping[read]=[{}, {}]
	if line.split("\t")[11]=="NM:i:0": mapping[read][0][contig_bins[contig]]=None
	if line.split("\t")[11]=="NM:i:0" or line.split("\t")[11]=="NM:i:1" or line.split("\t")[11]=="NM:i:2": mapping[read][1][contig_bins[contig]]=None


print "Going through forward and reverse reads and splitting them according to which bin they aligned to."
print "\nParsing forward reads and splitting them into files..."
ct=3
line1=""
line2=""
line3=""
line4=""
opened_files={}
for line in open(sys.argv[2]):
	ct+=1
	if ct==4:
		read=line[1:].split("/")[0].split()[0]
		ct=0

	if ct==0: line0=line
	elif ct==1: line1=line
	elif ct==2: line2=line
	elif ct==3:
		line3=line
		if read in mapping:
			for bin_name in mapping[read][0]:
				if bin_name+".strict" not in opened_files:
					opened_files[bin_name+".strict"] = {"ct": 0, "string": ""}
				opened_files[bin_name+".strict"]["string"]+=line0+line1+line2+line3
				opened_files[bin_name+".strict"]["ct"]+=1
				if opened_files[bin_name+".strict"]["ct"]==cache:
					f = open(sys.argv[4]+"/"+bin_name+".strict_1.fastq", 'a')
					f.write(opened_files[bin_name+".strict"]["string"])
					f.close()
					opened_files[bin_name+".strict"] = {"ct": 0, "string": ""}

			for bin_name in mapping[read][1]:
                                if bin_name+".permissive" not in opened_files:
                                        opened_files[bin_name+".permissive"] = {"ct": 0, "string": ""}
                                opened_files[bin_name+".permissive"]["string"]+=line0+line1+line2+line3
				opened_files[bin_name+".permissive"]["ct"]+=1
				if opened_files[bin_name+".permissive"]["ct"]==cache:
					f = open(sys.argv[4]+"/"+bin_name+".permissive_1.fastq", 'a')
					f.write(opened_files[bin_name+".permissive"]["string"])
					f.close()
					opened_files[bin_name+".permissive"] = {"ct": 0, "string": ""}

for bin_name in opened_files:
	f = open(sys.argv[4]+"/"+bin_name+"_1.fastq", 'a')
	f.write(opened_files[bin_name]["string"])
	f.close()


print "\nParsing reverse reads and splitting them into files..."
ct=3
line1=""
line2=""
line3=""
line4=""
opened_files={}
for line in open(sys.argv[3]):
        ct+=1
        if ct==4:
                read=line[1:].split("/")[0].split()[0]
                ct=0

        if ct==0: line0=line
        elif ct==1: line1=line
        elif ct==2: line2=line
        elif ct==3:
                line3=line
                if read in mapping:
	        	for bin_name in mapping[read][0]:
                                if bin_name+".strict" not in opened_files:
                                        opened_files[bin_name+".strict"] = {"ct": 0, "string": ""}
                                opened_files[bin_name+".strict"]["string"]+=line0+line1+line2+line3
                                opened_files[bin_name+".strict"]["ct"]+=1
                                if opened_files[bin_name+".strict"]["ct"]==cache:
                                        f = open(sys.argv[4]+"/"+bin_name+".strict_2.fastq", 'a')
                                        f.write(opened_files[bin_name+".strict"]["string"])
                                        f.close()
                                        opened_files[bin_name+".strict"] = {"ct": 0, "string": ""}

                        for bin_name in mapping[read][1]:
                                if bin_name+".permissive" not in opened_files:
                                        opened_files[bin_name+".permissive"] = {"ct": 0, "string": ""}
                                opened_files[bin_name+".permissive"]["string"]+=line0+line1+line2+line3
                                opened_files[bin_name+".permissive"]["ct"]+=1
                                if opened_files[bin_name+".permissive"]["ct"]==cache:
                                        f = open(sys.argv[4]+"/"+bin_name+".permissive_2.fastq", 'a')
                                        f.write(opened_files[bin_name+".permissive"]["string"])
                                        f.close()
                                        opened_files[bin_name+".permissive"] = {"ct": 0, "string": ""}

for bin_name in opened_files:
        f = open(sys.argv[4]+"/"+bin_name+"_2.fastq", 'a')
        f.write(opened_files[bin_name]["string"])
        f.close()

print "\nFinished splitting reads\n"
quit()

