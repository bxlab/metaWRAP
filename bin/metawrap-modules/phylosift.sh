#!/usr/bin/env bash

###########################################################################################################################################################
#       	                                                                                                                                          #
# This script is a wrapper aroung Phylosift - a taxonomic profiling software for reads.
#	                                                                                                                                 		  #
###########################################################################################################################################################

help_message () {
	echo ""
	echo "Usage: metaWRAP phylosift [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
	echo "NOTE: Phylosift does not work with some server scheduling services such as sbatch"
	echo "Options:"
	echo ""
	echo "	-1 STR          forward fastq reads"
	echo "	-2 STR          reverse fastq reads" 
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-s INT		number of reads to randomly sample (defualt=1000000)"
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
threads=1; out="false"; reads_1="false"; reads_2="false"
depth=1000000
# load in params
OPTS=`getopt -o ht:o:1:2:s: --long help -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
	case "$1" in
		-t) threads=$2; shift 2;;
		-o) out=$2; shift 2;;
		-1) reads_1=$2; shift 2;;
		-2) reads_2=$2; shift 2;;
		-s) depth=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
	esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ]; then 
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure the path to meta-scripts folder configured correctly in config.sh"
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

mkdir $out

# if sampling depth is specified, randomly subsample the fastq reads
comm "Subsampling down to $depth reads..."
paste $reads_1 $reads_2 | \
 awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | `#combine paired end reads onto one line` \
 shuf | head -n $depth | sed 's/\t\t/\n/g' | `#shuffle reads, select top N reads, and then restore tabulation` \
 awk '{print $1 > "'"${out}/subsample_1.fastq"'"; print $2 > "'"${out}/subsample_2.fastq"'"}' `#separate reads into F and R files`
reads_1=${out}/subsample_1.fastq
reads_2=${out}/subsample_2.fastq

announcement "Running Phylosift on $reads_1 and $reads_2..."

cmd="phylosift all --besthit --threads $threads --paired --output $out/phylosift_output --force $reads_1 $reads_2"
echo $cmd
$cmd

if [ ! -s $out/phylosift_output/subsample_1.fastq.xml ]; then 
	error "Phylosift did not finish correctly. Exiting."
else 
	cp $out/phylosift_output/subsample_1.fastq.html ${out}/kronagram.html
fi

########################################################################################################
########################              PHYLOSIFT PIPELINE COMPLETE!!!            ########################
########################################################################################################
announcement "PHYLOSIFT PIPELINE COMPLETE!!!"



