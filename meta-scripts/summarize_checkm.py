#!/usr/bin/env python
import sys
# This script summarizes the statistics of each bin by parsing 
# the checkm_folder/storage/bin_stats_ext.tsv file of the CheckM output


print "bin\tcompleteness\tcontamination\tGC\tlineage\tN50"
for line in open(sys.argv[1]):
	dic=eval(line.strip().split("\t")[1])

	#if dic["Completeness"]<20 or dic["Contamination"]>10: continue
	if "__" in dic["marker lineage"]: dic["marker lineage"]=dic["marker lineage"].split("__")[1]
	if "unbinned" in line.split("\t")[0]: name="Unbinned"
	else: name=line.split("\t")[0]


	print "\t".join([name, str(dic["Completeness"])[:5],\
	 str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
	 dic["marker lineage"], str(dic["N50 (contigs)"]) ])
