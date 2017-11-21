#!/bin/bash -l

##############################################################################################################################################################
#
# This script is meant to be a comprehensive solution for producing the best metagenomic assembly given paired end reads from one or more samples.
# Ideally it should take in fully QC'd reads. Due to its speed and efficiency, the default assembly will be done with Megahit, but if you have 
# sufficient resources (especially RAM) you may use metaSPAdes, which usually performs better. Short (<500bp) contigs are removed. The finall 
# assembly is then QCed by QUAST.
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
	echo "	-m INT          memory in GB (default=10)"
	echo "	-t INT          number of threads (defualt=1)"
	echo ""
	echo "	--metaspades	assemble with metaspades instead of megahit"
	echo "";}
comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }



########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
source config-metawrap

# default params
mem=10; threads=1; out="false"; reads_1="false"; reads_2="false"
# long options defaults
metaspades_assemble=false


# load in params
OPTS=`getopt -o ht:m:o:1:2: --long help,metaspades -- "$@"`
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
		-h | --help) help_message; exit 1; shift 1;;
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


mkdir $out

if [ "$metaspades_assemble" = true ]; then
	########################################################################################################
	########################               ASSEMBLING WITH METASPADES               ########################
	########################################################################################################
        announcement "ASSEMBLING WITH METASPADES"
	
	comm "Using reads $reads_1 and $reads_2 for assembly,"
	if [ ! -f "$reads_1" ] || [ ! -f "$reads_2" ]; then error "Read files $reads_1 and/or $reads_2 dont exist. Exiting."; fi
	
	metaspades.py --tmp-dir ${out}/metaspades.tmp -t $threads -m $mem -o ${out}/metaspades -1 $reads_1 -2 $reads_2 
	rm -r ${out}/metaspades.tmp
	if [ ! -f "${out}/metaspades/scaffolds.fasta" ]; then error "Something went wrong with metaSPAdes assembly. Exiting."; fi

else

	########################################################################################################
	########################                ASSEMBLING WITH MEGAHIT                 ########################
	########################################################################################################
        announcement "ASSEMBLING READS WITH MEGAHIT"
	
	rm -r ${out}/megahit
	comm "assembling $reads_1 and $reads_2 with megahit"
	megahit -1 $reads_1 -2 $reads_2 -o ${out}/megahit -t $threads -m ${mem}000000000

	if [ ! -f "${out}/megahit/final.contigs.fa" ]; then error "Something went wrong with reassembling with Megahit. Exiting."; fi

fi

########################################################################################################
########################         COMBINE AND FORMAT THE TWO ASSEMBLIES          ########################
########################################################################################################
announcement "COMBINE AND FORMAT THE TWO ASSEMBLIES"

if [ "$metaspades_assemble" = true ]; then
	comm "Reads were assembled with metaspades only"
	${SOFT}/rm_short_contigs.py 500 ${out}/metaspades/scaffolds.fasta > ${out}/final_assembly.fasta
	if [ ! -s ${out}/final_assembly.fasta ]; then error "metaspades failed to produce long contigs (>500bp). Exiting."; fi
elif [ "$metaspades_assemble" = false ]; then 
	comm "Reads were assembled with megahit only"
	${SOFT}/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa > ${out}/megahit/fixed.contigs.fa
	${SOFT}/rm_short_contigs.py 500 ${out}/megahit/fixed.contigs.fa > ${out}/megahit/long.contigs.fa
	if [ ! -s ${out}/megahit/long.contigs.fa ]; then error "megahit failed to produce long (>500bp) contigs. Exiting."; fi
	${SOFT}/sort_contigs.py ${out}/megahit/long.contigs.fa > ${out}/final_assembly.fasta
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

