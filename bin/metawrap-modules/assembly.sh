#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script is meant to be a comprehensive solution for producing the best metagenomic assembly given paired end reads from one or more samples.
# Ideally it should take in fully QC'd reads. First, the reads are assembled with metaSPAdes3.10, then all the reads that did not map back to the
# contigs are re-assembled with MEGAHIT (which works better on lower coverage contigs. The resulting assemblies are combined, sorted, and short 
# contigs are removed. The finall assembly is then QCed by QUAST.
#
# Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################



help_message () {
	echo ""
	echo "Usage: metaWRAP assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
	echo "Options:"
	echo ""
	echo "	-1 STR          forward fastq reads"
	echo "	-2 STR          reverse fastq reads" 
	echo "	-o STR          output directory"
	echo "	-m INT          memory in GB (default=24)"
	echo "	-t INT          number of threads (defualt=1)"
	echo "	-l INT		minimum length of assembled contigs (default=1000)"
	echo ""
	echo "	--megahit	assemble with megahit (default)"
	echo "	--metaspades	assemble with metaspades instead of megahit (better results but slower and higher memory requirement)"
	echo "";}
comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }



########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
config_file=$(which config-metawrap)
source $config_file

# default params
mem=24; threads=1; out="false"; reads_1="false"; reads_2="false"; min_len=1000
# long options defaults
metaspades_assemble=false; megahit_assemble=true


# load in params
OPTS=`getopt -o ht:m:o:1:2:l: --long help,metaspades,megahit -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
	case "$1" in
		-t) threads=$2; shift 2;;
		-m) mem=$2; shift 2;;
		-o) out=$2; shift 2;;
		-1) reads_1=$2; shift 2;;
		-2) reads_2=$2; shift 2;;
		-l) min_len=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		--megahit) megahit_assemble=true; shift 1;;
		--metaspades) metaspades_assemble=true; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ] || [ "$threads" = "false" ] || [ "$mem" = "false" ]; then 
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure the meta-scripts path is configured correctly in config.sh file"
fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################


if [ ! -d $out ]; then
        mkdir $out;
else
        echo "Warning: $out already exists."
fi

if [ "$metaspades_assemble" = true ]; then
	########################################################################################################
	########################               ASSEMBLING WITH METASPADES               ########################
	########################################################################################################
        announcement "ASSEMBLING WITH METASPADES"
	
	comm "Using reads $reads_1 and $reads_2 for assembly,"
	if [ ! -f "$reads_1" ] || [ ! -f "$reads_2" ]; then error "Read files $reads_1 and/or $reads_2 dont exist. Exiting."; fi
	
	mkdir ${out}/metaspades.tmp

	if [[ -s ${out}/metaspades/spades.log ]]; then
		metaspades.py -o ${out}/metaspades --restart-from last -t $threads -m $mem --tmp-dir ${out}/metaspades.tmp
	else
		metaspades.py --tmp-dir ${out}/metaspades.tmp -t $threads -m $mem -o ${out}/metaspades -1 $reads_1 -2 $reads_2
	fi

	if [ ! -f "${out}/metaspades/scaffolds.fasta" ]; then error "Something went wrong with metaSPAdes assembly. Exiting."; fi
	rm -r ${out}/metaspades.tmp
fi


if [ "$megahit_assemble" = true ] && [ "$metaspades_assemble" = true ]; then
	########################################################################################################
	########################              SORTING OUT UNASSEMBLED READS             ########################
	########################################################################################################
        announcement "SORTING OUT UNASSEMBLED READS"

	#take only the scaffolds over 1.5kb (the minimum needed for binning with metaBAT):
	${SOFT}/rm_short_contigs.py 1500 ${out}/metaspades/scaffolds.fasta > ${out}/metaspades/long_scaffolds.fasta
	if [ ! -s ${out}/metaspades/long_scaffolds.fasta ]; then
		error "metaSPAdes produced no contigs over 1.5kb, so the rest of the pipeline wont work. Exiting."; fi


	bwa index ${out}/metaspades/long_scaffolds.fasta
	
	#sort out and store reads that dont map back to the assembly:
	bwa mem -t $threads ${out}/metaspades/long_scaffolds.fasta $reads_1 $reads_2 | grep -v NM:i: > ${out}/unused_by_metaspades.sam
	${SOFT}/sam_to_fastq.py ${out}/unused_by_metaspades.sam > ${out}/unused_by_metaspades.fastq
	rm ${out}/unused_by_metaspades.sam
	
	if [[ ! -s ${out}/unused_by_metaspades.fastq ]]; then error "Something went wrong with pulling out unassembled reads. Exiting."; fi
fi

if [ "$megahit_assemble" = true ]; then
	########################################################################################################
	########################        ASSEMBLING REMAINDER READS WITH MEGAHIT         ########################
	########################################################################################################
        announcement "ASSEMBLING READS WITH MEGAHIT"
	
	#rm -r ${out}/megahit
	mkdir ${out}/megahit.tmp
	if [ "$metaspades_assemble" = true ]; then
		comm "assembling ${out}/unused_by_metaspades.fastq with megahit"
		megahit\
		 -r ${out}/unused_by_metaspades.fastq\
		 -o ${out}/megahit\
		 -t $threads\
		 -m ${mem}000000000\
		 --tmp-dir ${out}/megahit.tmp
		mv ${out}/unused_by_metaspades.fastq ${out}/metaspades/
	else
		comm "assembling $reads_1 and $reads_2 with megahit"
		megahit\
		 -1 $reads_1 -2 $reads_2\
		 -o ${out}/megahit\
		 --tmp-dir ${out}/megahit.tmp\
		 -t $threads -m ${mem}000000000\
		 --continue
	fi

	if [ ! -f "${out}/megahit/final.contigs.fa" ]; then error "Something went wrong with reassembling with Megahit. Exiting."; fi
	rm -r ${out}/megahit.tmp

fi

########################################################################################################
########################         COMBINE AND FORMAT THE TWO ASSEMBLIES          ########################
########################################################################################################
announcement "FORMAT THE ASSEMBLY"

if [ "$megahit_assemble" = true ] && [ "$metaspades_assemble" = true ]; then
	comm "Reads were assembled with metaspades AND megahit"
	${SOFT}/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa $min_len > ${out}/megahit/long.contigs.fa
	
	cp ${out}/metaspades/long_scaffolds.fasta ${out}/combined_assembly.fasta
	cat ${out}/megahit/long.contigs.fa >> ${out}/combined_assembly.fasta
	${SOFT}/sort_contigs.py ${out}/combined_assembly.fasta > ${out}/final_assembly.fasta
	rm ${out}/combined_assembly.fasta
elif [ "$metaspades_assemble" = true ]; then
	comm "Reads were assembled with metaspades"
	${SOFT}/rm_short_contigs.py $min_len ${out}/metaspades/scaffolds.fasta > ${out}/final_assembly.fasta
	if [ ! -s ${out}/final_assembly.fasta ]; then error "metaspades failed to produce long contigs (>500bp). Exiting."; fi
elif [ "$megahit_assemble" = true ]; then 
	comm "Reads were assembled with megahit. Formatting and sorting the assembly..."
	${SOFT}/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa $min_len > ${out}/final_assembly.fasta
	if [[ $? -ne 0 ]]; then error "something went wrong with formating and sorting the assembly. Exiting..."; fi
	if [ ! -s ${out}/final_assembly.fasta ]; then error "megahit failed to produce long (>${min_len}) contigs. Exiting."; fi
else
	error "Reads were NOT assembled with metaspades OR megahit..."
fi

if [[ ! -s ${out}/final_assembly.fasta ]]; then error "Something went wrong with joining the assemblies. Exiting."; fi

########################################################################################################
########################             RUNNING ASSEMBLY QC WITH QUAST             ########################
########################################################################################################
announcement "RUNNING ASSEMBLY QC WITH QUAST"

quast -t $threads -o ${out}/QUAST_out -m 500 ${out}/final_assembly.fasta
cp ${out}/QUAST_out/report.html ${out}/assembly_report.html

if [[ ! -s ${out}/assembly_report.html ]]; then error "Something went wrong with running QUAST. Exiting."; fi

########################################################################################################
########################      ASSEMBLY PIPELINE COMPLETED SUCCESSFULLY!!!       ########################
########################################################################################################
announcement "ASSEMBLY PIPELINE COMPLETED SUCCESSFULLY!!!"

