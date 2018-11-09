#!/usr/bin/env python2.7
import sys

# This script takes in the reassembled_bins.stats file of the binning module and choses the best possible
# case scenerio for each bin. Each bin should have 3 versions of itself: the original bin (.orig), the
# "strict" reassembled bin (.strict), and the "permissive" reassembled bin (permissive). This program
# decides which one of the three is better for each bin. The goal is to have the highest completion, while
# maintaining low contamination. To do so, we compute a "score" for each assembly. SCORE=Compl + (100-Cont)*5
# Note that we put more value into contamination, as its usually better to have lower contamination than
# slightly higher completion. Ambiguities are settled by the N50 score.

best_bins={}
for line in open(sys.argv[1]):
	if "contamination" in line: continue
	cut=line.strip().split("\t")

	# dont consider bins with < min completion and > max contamination
	if float(cut[1])<float(sys.argv[2]) or float(cut[2])>float(sys.argv[3]): continue

	bin_name=".".join(cut[0].split(".")[:-1])	
	style=cut[0].split(".")[-1]
	score = float(cut[1]) + 5*(100-float(cut[2]))
	n50 = int(cut[5])
	if bin_name not in best_bins: 
		best_bins[bin_name]=[style, score, n50]
	else:
		if score>best_bins[bin_name][1]:
			best_bins[bin_name]=[style, score, n50]
		if score==best_bins[bin_name][1] and n50>best_bins[bin_name][2]:
			best_bins[bin_name]=[style, score, n50]
		if score<best_bins[bin_name][1]:
			continue
	


for i in best_bins:
	print i+'.'+best_bins[i][0]
	


