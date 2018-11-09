#!/usr/bin/env python2.7
import sys, os

'''
This script takes in two folders containing bins from two binning methods, and their complementary 
contamination and completion scores derived from CheckM, then matches corresponding bins in the two
sets based on a minimum of 80% overlap (by length) of the bins, and finally, decides which of the two
bin versions is best in that particular bin.

Then it makes a new folder into which it puts the best version of each bin (changing the naming in
the process), and also makes a new .stats file which is consistant with the new bin folder.

Usage:
./script bin_folder_1 bin_folder_2 stats_file_1 stats_file_2 output_folder min_completion max_contaminaiton
'''

c=float(sys.argv[6])
x=float(sys.argv[7])

# load a list of good bins (>70% complete, <10% contaminated) to save time (wont look at bad bins later on).
print "Loading list of good bins (comp>" + str(c)+ "%, cont<" + str(x) + "%)"
good_bins_1={}
good_bins_2={}
for line in open(sys.argv[3]):
	if "completeness" in line: continue
	cut=line.strip().split('\t')
	if float(cut[1])>c and float(cut[2])<x: good_bins_1[cut[0]+'.fa']=None
for line in open(sys.argv[4]):
	if "completeness" in line: continue
        cut=line.strip().split('\t')
        if float(cut[1])>c and float(cut[2])<x: good_bins_2[cut[0]+'.fa']=None


# these are the dictionaries storing the contig names and lengths of each bin
bins_1={}
bins_2={}

'''
print  "load in the info about the contigs in each bin..."
for bin_file in good_bins_1:
	bins_1[bin_file]={}
	for line in open(sys.argv[1] + '/' + bin_file):
		if not line.startswith('>'): continue
		bins_1[bin_file][line[1:-1]] = int(line.split('_')[3])
for bin_file in good_bins_2:
        bins_2[bin_file]={}
        for line in open(sys.argv[2] + '/' + bin_file):
                if not line.startswith('>'): continue
                bins_2[bin_file][line[1:-1]] = int(line.split('_')[3])	
'''


print  "load in the info about the contigs in each bin..."
for bin_file in good_bins_1:
        bins_1[bin_file]={}
	contig_len=0
	contig_name=""
        for line in open(sys.argv[1] + '/' + bin_file):
                if not line.startswith('>'):
			contig_len+=len(line.strip())
                else:
			if contig_name!="":
				bins_1[bin_file][contig_name] = contig_len
				contig_len=0
			contig_name=line[1:-1]
	bins_1[bin_file][contig_name] = contig_len

for bin_file in good_bins_2:
	contig_len=0
	contig_name=""
        bins_2[bin_file]={}
        for line in open(sys.argv[2] + '/' + bin_file):
                if not line.startswith('>'): 
			contig_len+=len(line.strip())
                else:
			if contig_name!="":
				bins_2[bin_file][contig_name] = contig_len
				contig_len=0
			contig_name=line[1:-1]
	bins_2[bin_file][contig_name] = contig_len




print "make all bossible comparisons between the two bin sets, and record total % idential length"
all_bin_pairs={}
for bin_1 in good_bins_1:
	all_bin_pairs[bin_1]={}
	for bin_2 in good_bins_2:
		# find idential contigs between bin_1 and bin_2
		match_1_length=0
		match_2_length=0
		mismatch_1_length=0
		mismatch_2_length=0

		for contig in bins_1[bin_1]:
			if contig in bins_2[bin_2]: match_1_length+=bins_2[bin_2][contig]
			else: mismatch_1_length+=bins_1[bin_1][contig]
		for contig in bins_2[bin_2]:
                        if contig in bins_1[bin_1]: match_2_length+=bins_1[bin_1][contig]
			else: mismatch_2_length+=bins_2[bin_2][contig]

		# chose the highest % ID, dependinsh of which bin is  asubset of the other
		ratio_1=100*match_1_length/(match_1_length+mismatch_1_length)
		ratio_2=100*match_2_length/(match_2_length+mismatch_2_length)

		all_bin_pairs[bin_1][bin_2]=max([ratio_1, ratio_2])


print "load in completion and contamination scores of all the bins"
bins_1_stats={}
bins_2_stats={}
bins_1_summary={}
bins_2_summary={}
for line in open(sys.argv[3]):
	if "completeness" in line: 
		bins_1_summary["header"]=line
		continue
	cut=line.strip().split('\t')
	bins_1_stats[cut[0]+'.fa']=(float(cut[1]), float(cut[2]))
	bins_1_summary[cut[0]+'.fa']=line

for line in open(sys.argv[4]):
        if "completeness" in line: 
		bins_2_summary["header"]=line
		continue
        cut=line.strip().split('\t')
        bins_2_stats[cut[0]+'.fa']=(float(cut[1]), float(cut[2]))
	bins_2_summary[cut[0]+'.fa']=line



# go through all good bins and chose best ones
print "go through first group, pull out identical bins from second group, and choose best"
os.system("mkdir "+sys.argv[5])
new_summary_file=bins_1_summary["header"]
bins_2_matches={}
bin_ct=1
for bin_1 in all_bin_pairs:
	score=bins_1_stats[bin_1][0] - bins_1_stats[bin_1][1]*5
	found_better=False
	for bin_2 in all_bin_pairs[bin_1]:
		# check for sufficient overlap (80% bin length)
		if all_bin_pairs[bin_1][bin_2]<80: continue
		bins_2_matches[bin_2]=None
		# check if this bin is better than original
		if (bins_2_stats[bin_2][0]-bins_2_stats[bin_2][1]*5) > score:
			cmd = "cp " + sys.argv[2] + '/' + bin_2 + " " + sys.argv[5] + "/bin." + str(bin_ct) + ".fa"
			new_summary_file+="bin." + str(bin_ct) + "\t" + "\t".join(bins_2_summary[bin_2].split("\t")[1:])
			found_better=True
	if found_better==False: 
		new_summary_file+="bin." + str(bin_ct) + "\t" + "\t".join(bins_1_summary[bin_1].split("\t")[1:])
		cmd = "cp " + sys.argv[1] + '/' + bin_1 + " " + sys.argv[5] + "/bin." + str(bin_ct) + ".fa"
	os.system(cmd)
	bin_ct+=1

print "retrieve bins from second group that were not found in first group"
for bin_2 in bins_2_stats:
	if bins_2_stats[bin_2][0]<c or bins_2_stats[bin_2][1]>x: continue
	if bin_2 in bins_2_matches: continue
	new_summary_file+="bin." + str(bin_ct) + "\t" + "\t".join(bins_2_summary[bin_2].split("\t")[1:])
	cmd = "cp " + sys.argv[2] + '/' + bin_2 + " " + sys.argv[5] + "/bin." + str(bin_ct) + ".fa"
	os.system(cmd)
	bin_ct+=1


f = open(sys.argv[5]+".stats", 'w') 
f.write(new_summary_file) 

print "There were " + str(bin_ct) + " bins cherry-picked from the original sets!"


