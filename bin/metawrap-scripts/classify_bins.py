#!/usr/bin/env python2.7
import sys,os

def add_to_tree( tree, tax_list, length ):
	if len(tax_list)==0: return tree 
	else:
		if tax_list[0] not in tree: tree[tax_list[0]]=[length, {}]
		else: tree[tax_list[0]][0]+=length
		add_to_tree( tree[tax_list[0]][1], tax_list[1:], length )
	return tree

def traverse(tree, taxonomy, weight):
	if len(tree)==0:
		return taxonomy
	total_score=0
	max_score=0
	max_class=""
	for k in tree:
		total_score+=tree[k][0]
		if tree[k][0]>max_score: 
			max_score=tree[k][0]
			max_class=k
	if weight!=0: total_score=weight

	if 100*max_score/total_score>50:
		taxonomy.append(max_class) #+"_"+str(100*max_score/total_score))
		taxonomy=traverse(tree[max_class][1], taxonomy, tree[max_class][0])
	else:
		return taxonomy
	return taxonomy


# load in classifications of each contig
taxonomy={}
for line in open (sys.argv[1]):
	cut=line.strip().split("\t")
	if len(cut)<2: continue
	taxonomy[cut[0]] = cut[1]


# loop over bins
for filename in os.listdir(sys.argv[2]):
	tax_tree={}
	length=0
	for line in open(sys.argv[2]+'/'+filename):
		if line[0]=='>' and length==0:
			contig=line[1:-1]
		elif line[0]=='>' and length>0:
			if contig in taxonomy: 
				#print "\t".join([contig, str(length), taxonomy[contig]])
				taxcut=taxonomy[contig].split(';')
				tax_tree = add_to_tree(tax_tree, taxcut, length)
			contig=line[1:-1]
			length=0
		else:
			length+=len(line.strip())
	if contig in taxonomy:
		taxcut=taxonomy[contig].split(';')
		tax_tree = add_to_tree(tax_tree, taxcut, length)

	#print tax_tree
	consensus=traverse(tax_tree, [], 0)
	print filename + "\t" + ";".join(consensus)

