#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script takes in a set of bins and any number of paired-end read sets from metagenomic samples, and estimates the abundance of each bin in each 
# sample with salmon. It then uses Seaborn to make clustered heatmaps genome abundances.
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: metaWRAP quant_bins [options] -b bins_folder -o output_dir -a assembly.fa readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]"
	echo "Options:"
	echo ""
	echo "	-b STR          folder containing draft genomes (bins) in fasta format"
	echo "	-o STR          output directory"
	echo "	-a STR		fasta file with entire metagenomic assembly (strongly recommended!)"
	echo "	-t INT		number of threads"
	echo ""
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
threads=1; out=false; bin_folder=false; assembly=false
# long options defaults

# load in params
OPTS=`getopt -o ht:o:b:a: --long help -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
                -o) out=$2; shift 2;;
                -b) bin_folder=$2; shift 2;;
		-a) assembly=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ $out = false ] || [ $bin_folder = false ]; then 
	comm "Non-optional parameters -b and -o were not entered"
	help_message; exit 1
fi

# check for at least one pair of read fastq files:
F="no"; R="no"
for num in "$@"; do
	if [[ $num == *"_1.fastq" ]]; then F="yes"; fi
	if [[ $num == *"_2.fastq" ]]; then R="yes"; fi
done
if [ $F = "no" ] || [ $R = "no" ]; then
	comm "Unable to find proper fastq read pair in the format *_1.fastq and *_2.fastq"
	help_message; exit 1
fi

#determine number of fastq read files provided:
num_of_F_read_files=$(for I in "$@"; do echo $I | grep _1.fastq; done | wc -l)
num_of_R_read_files=$(for I in "$@"; do echo $I | grep _2.fastq; done | wc -l)

comm "$num_of_F_read_files forward and $num_of_R_read_files reverse read files detected"
if [ ! $num_of_F_read_files == $num_of_R_read_files ]; then error "Number of F and R reads must be the same!"; fi


# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################


########################################################################################################
########################         SETTING UP OUTPUT AND INDEXING ASSEMBLY        ########################
########################################################################################################
announcement "SETTING UP OUTPUT AND INDEXING ASSEMBLY"

# setting up the output folder
if [ ! -d ${out}/all_bin_contigs ]; then
        mkdir $out
else
        warning "Warning: $out already exists."
fi

if [ $assembly = false ]; then
	comm "Concatinating bins into a metagenomic assembly file."
	if [[ -s ${out}/assembly.fa ]]; then rm ${out}/assembly.fa; fi
	cat $bin_folder/* >> ${out}/assembly.fa
	assembly=${out}/assembly.fa
else
	if [ ! -f $assembly ]; then error "Assembly file $assembly does not exist. Exiting..."; fi
fi

# Index the assembly
comm "Indexing assembly file with salmon. Ignore any warnings"
salmon index -p $threads -t $assembly -i ${out}/assembly_index
if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi



########################################################################################################
########################      ALIGNING READS FROM ALL SAMPLES BACK TO BINS      ########################
########################################################################################################
announcement "ALIGNING READS FROM ALL SAMPLES BACK TO BINS WITH SALMON"

mkdir ${out}/alignment_files
# If there are several pairs of reads passed, they are processed sepperately
for arg in "$@"; do
	if [[ $arg == *"_1.fastq"* ]]; then
		reads_1=$arg
		reads_2=${arg%_*}_2.fastq
		if [ "$reads_2" = "false" ]; then error "$reads_2 does not exist. Exiting..."; fi

		tmp=${reads_1##*/}
		sample=${tmp%_*}

		comm "processing sample $sample with reads $reads_1 and $reads_2..."
		salmon quant -i ${out}/assembly_index --libType IU -1 $reads_1 -2 $reads_2 -o ${out}/alignment_files/${sample}.quant --meta -p $threads
		if [[ $? -ne 0 ]]; then error "Something went wrong with aligning $sample fastq files back to assembly! Exiting..."; fi
	fi
done

comm "summarize salmon files..."
home=$(pwd)
cd ${out}/alignment_files/
${SOFT}/summarize_salmon_files.py
if [[ $? -ne 0 ]]; then error "something went wrong with summarizing salmon output with python script! Exiting..."; fi
cd $home
mkdir ${out}/quant_files
for f in $(ls ${out}/alignment_files/ | grep .quant.counts); do mv ${out}/alignment_files/$f ${out}/quant_files/; done



########################################################################################################
########################        EXTRACTING AVERAGE ABUNDANCE OF EACH BIN        ########################
########################################################################################################
announcement "EXTRACTING AVERAGE ABUNDANCE OF EACH BIN"

n=$(ls ${out}/quant_files/ | grep counts | wc -l)
if [[ $n -lt 1 ]]; then error "There were no files found in ${out}/quant_files/"; fi
comm "There were $n samples detected. Making abundance table!"

${SOFT}/split_salmon_out_into_bins.py ${out}/quant_files/ $bin_folder $assembly > ${out}/bin_abundance_table.tab
if [[ $? -ne 0 ]]; then error "something went wrong with making summary abundance table. Exiting..."; fi
comm "Average bin abundance table stored in ${out}/abundance_table.tab"




########################################################################################################
########################            MAKING GENOME ABUNDANCE HEATMAP             ########################
########################################################################################################
if [[ $n -gt 1 ]]; then
	announcement "MAKING GENOME ABUNDANCE HEATMAP WITH SEABORN"

	comm "making heatmap with Seaborn"
	${SOFT}/make_heatmap.py ${out}/bin_abundance_table.tab ${out}/bin_abundance_heatmap.png
	if [[ $? -ne 0 ]]; then error "something went wrong with making the heatmap. Exiting..."; fi

	comm "cleaning up..."
	rm -r ${out}/alignment_files/ 
else
	warning "Cannot make clustered heatmap with just one sample... Skipping heatmap"
fi

########################################################################################################
########################     QUANT_BINS PIPELINE SUCCESSFULLY FINISHED!!!       ########################
########################################################################################################
announcement "QUANT_BINS PIPELINE SUCCESSFULLY FINISHED!!!"

