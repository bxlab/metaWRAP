#!/usr/bin/env bash

##############################################################################################################################################################
#
# Runs CheckM on a directory full of MAGs/bins
#
##############################################################################################################################################################


source $(which config-metawrap)

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }

run_checkm () {
	comm "Running CheckM on ${1} bins"

	p_threads=40
	threads=40

	mkdir ${1}.tmp
	checkm lineage_wf -x fa ${1} ${1}.checkm -t $threads --tmpdir ${1}.tmp --pplacer_threads $p_threads
	if [[ ! -s ${1}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r ${1}.tmp
	${SOFT}/summarize_checkm.py ${1}.checkm/storage/bin_stats_ext.tsv ${1##*/}\
	| (read -r; printf "%s\n" "$REPLY"; sort -rn -k2) > ${1}.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
}


run_checkm $1

