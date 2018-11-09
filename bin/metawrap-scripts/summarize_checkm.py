#!/usr/bin/env python2.7
import sys
# This script summarizes the statistics of each bin by parsing 
# the checkm_folder/storage/bin_stats_ext.tsv file of the CheckM output


if len(sys.argv)==3: 
	binner=sys.argv[2]
	print "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner"
elif len(sys.argv)==4:
	source={}
	for line in open(sys.argv[3]):
		cut=line.strip().split("\t")
		source[cut[0]]=cut[7]
	print "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\tbinner"
else:
	print "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize"


for line in open(sys.argv[1]):
	dic=eval(line.strip().split("\t")[1])

	#if dic["Completeness"]<20 or dic["Contamination"]>10: continue
	if "__" in dic["marker lineage"]: dic["marker lineage"]=dic["marker lineage"].split("__")[1]
	#if "unbinned" in line.split("\t")[0]: name="Unbinned"
	name=line.split("\t")[0]


	if len(sys.argv)==3:	
		print "\t".join([name, str(dic["Completeness"])[:5],\
		 str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
		 dic["marker lineage"], str(dic["N50 (contigs)"]),\
		 str(dic["Genome size"]), binner])

	elif len(sys.argv)==4:
		print "\t".join([name, str(dic["Completeness"])[:5],\
		 str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
		 dic["marker lineage"], str(dic["N50 (contigs)"]),\
		 str(dic["Genome size"]), source[name]])

	else:
		print "\t".join([name, str(dic["Completeness"])[:5],\
		 str(dic["Contamination"])[:5], str(dic["GC"])[:5],\
		 dic["marker lineage"], str(dic["N50 (contigs)"]),\
		 str(dic["Genome size"])])
