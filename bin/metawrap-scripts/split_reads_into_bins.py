#!/usr/bin/env python2.7
# Usage: ./script mapping_file.dict reads_1.fastq reads_2.fastq output_dir
import sys, os

print "\nLoading in dictionary containing bin mapping of all the reads...\n"
mapping=eval(open(sys.argv[1]).readline())


print "\n\nParsing forward reads and splitting them into files..."
ct=3
line1=""
line2=""
line3=""
line4=""
opened_files={}
for line in open(sys.argv[2]):
	ct+=1
	if ct==4:
		read=line[1:].split("/")[0]
		ct=0
	
	if ct==0: line0=line
	elif ct==1: line1=line
	elif ct==2: line2=line
	elif ct==3: 
		line3=line
		if read in mapping:
			for bin_name in mapping[read][0]:
				if bin_name+"_strict" not in opened_files: 
					print "Opening file "+sys.argv[4]+'/'+bin_name+".strict_1.fastq"
					opened_files[bin_name+"_strict"] = open(sys.argv[4]+"/"+bin_name+".strict_1.fastq", 'w')
				opened_files[bin_name+"_strict"].write(line0)
				opened_files[bin_name+"_strict"].write(line1)
				opened_files[bin_name+"_strict"].write(line2)
				opened_files[bin_name+"_strict"].write(line3)
                        
			for bin_name in mapping[read][1]:
                                if bin_name+"permissive" not in opened_files:
                                        print "Opening file "+sys.argv[4]+'/'+bin_name+".permissive_1.fastq"
                                        opened_files[bin_name+"permissive"] = open(sys.argv[4]+"/"+bin_name+".permissive_1.fastq", 'w')
                                opened_files[bin_name+"permissive"].write(line0)
                                opened_files[bin_name+"permissive"].write(line1)
                                opened_files[bin_name+"permissive"].write(line2)
                                opened_files[bin_name+"permissive"].write(line3)
print "Closing bin fastq files"
for f in opened_files: opened_files[f].close()


print "\n\nParsing reverse reads and splitting them into files..."
ct=3
line1=""
line2=""
line3=""
line4=""
opened_files={}
for line in open(sys.argv[3]):
        ct+=1
        if ct==4:
                read=line[1:].split("/")[0]
                ct=0

        if ct==0: line0=line
        elif ct==1: line1=line
        elif ct==2: line2=line
        elif ct==3:
                line3=line
                if read in mapping:
                        for bin_name in mapping[read][0]:
                                if bin_name+"_strict" not in opened_files:
                                        print "Opening file "+sys.argv[4]+'/'+bin_name+".strict_2.fastq"
                                        opened_files[bin_name+"_strict"] = open(sys.argv[4]+"/"+bin_name+".strict_2.fastq", 'w')
                                opened_files[bin_name+"_strict"].write(line0)
                                opened_files[bin_name+"_strict"].write(line1)
                                opened_files[bin_name+"_strict"].write(line2)
                                opened_files[bin_name+"_strict"].write(line3)

                        for bin_name in mapping[read][1]:
                                if bin_name+"permissive" not in opened_files:
                                        print "Opening file "+sys.argv[4]+'/'+bin_name+".permissive_2.fastq"
                                        opened_files[bin_name+"permissive"] = open(sys.argv[4]+"/"+bin_name+".permissive_2.fastq", 'w')
                                opened_files[bin_name+"permissive"].write(line0)
                                opened_files[bin_name+"permissive"].write(line1)
                                opened_files[bin_name+"permissive"].write(line2)
                                opened_files[bin_name+"permissive"].write(line3)

print "Closing bin fastq files"
for f in opened_files: opened_files[f].close()
				
				
