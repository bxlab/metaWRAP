#!/bin/bash -l

##############################################################################################################################################################
#
# This script is meant to be run on paired end reads (with extensions *_1.fastq and *_2.fastq) or assembled contigs (*.fa or *.fasta).
# The script runs KRAKEN on the sequences, then translates them to taxonomy form with kraken-translate. Then in-house scripts are used to 
# parse out the taxonomy into the format for KRONA-TOOLS, colapse the file to save memory, and finally produce a prety kronagram with all the files.
#
# NOTE: KRAKEN and KronaTools requires instalation, and be sure the configure the right path to the KRAKEN folder and the KRAKEN database.
#
# Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################


help_message () {
	echo "Usage: run on any number of fasta assembly files and/or or paired-end reads."
	echo "./kraken.sh [options] -o output_dir assembly.fasta reads_1.fastq reads_2.fastq ..."
	echo "Options:"
	echo "" 
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads"
	echo "	-d INT		read sampling depth (default=all)"
	echo "";}

# function to print out error messages
error () {
	echo "************************************************************************************************************************"
	echo "*****************************************            ERROR!           **************************************************"
	echo " "; echo $1; echo " "
	echo "************************************************************************************************************************"; exit 1; }
# function to print out warning messages
warning () {
	echo "************************************************************************************************************************"
	echo "*****************************************           WARNING!          **************************************************"
	echo " "; echo $1; echo " "
	echo "************************************************************************************************************************";}
# funciton to print out comments throught the pipeline
comm () {
	echo ""; echo "----------------------------------------------------------------------------------------------"
	echo $1
	echo "----------------------------------------------------------------------------------------------"; echo ""; }



########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
source ${0%/*}/config.sh


# Set defaults
threads=1; out="false"; depth="all"
# Load in options
while getopts ht:o:d: option; do
	case "${option}" in
		h) help_message; exit 1;;
		t) threads=${OPTARG};;
		o) out=${OPTARG};;
		d) depth=${OPTARG};;
	esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$#" -lt 1 ] ; then 
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure the meta-scriptss folder is in the same folder as the other reads_qc"
fi

# Checks for KRAKEN database 
if [ ! -d "$KRAKEN_DB" ]; then
	error "The folder $KRAKEN_DB doesnt exist. Please look and the the script and manually set the correct path to the KRAKEN standard database. \
	If you do not have it yet, you can download it like this: ${KRAKEN}/kraken-build --standard --threads $threads --db ${KRAKEN}/kraken/KRAKEN_DB"
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################




echo " "
echo "########################################################################################################"
echo "########################              RUNNING KRAKEN ON ALL FILES               ########################"
echo "########################################################################################################"
echo " "

# setting up the output folder
mkdir $out

# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	#process fastq files
	if [[ $num == *"_1.fastq" ]]; then
		reads_1=$num
		reads_2=${num%_*}_2.fastq
		if [ "$reads_2" = "false" ]; then error "$reads_2 does not exist. Exiting..."; fi
		
		tmp=${reads_1##*/}
		sample=${tmp%_*}
		comm "Now processing $reads_1 and $reads_2"

		# if sampling depth is specified, randomly subsample the fastq reads
		if [ ! "$depth" = "all" ]; then
			comm "subsampling $depth reads" 
			paste $reads_1 $reads_2 | \
			 awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | `#combine paired end reads onto one line` \
			 shuf | head -n $depth | sed 's/\t\t/\n/g' | `#shuffle reads, select top N reads, and then restore tabulation` \
			 awk '{print $1 > "'"${out}/tmp_1.fastq"'"; print $2 > "'"${out}/tmp_2.fastq"'"}' `#separate reads into F and R files`
			reads_1=${out}/tmp_1.fastq
			reads_2=${out}/tmp_2.fastq
		fi

		#run kraken	
		kraken --db ${KRAKEN_DB} --fastq-input --paired --threads $threads\
	 	--output ${out}/${sample}.krak $reads_1 $reads_2

		if [ ! "$depth" = "all" ]; then rm ${out}/tmp_1.fastq ${out}/tmp_2.fastq; fi
	fi

	#process fasta files
	if [[ $num == *"fa" ]] || [[ $num == *"fasta" ]]; then
		tmp=${num##*/}
		sample=${tmp%.*}
		comm "Now processing $num"

		kraken --db ${KRAKEN_DB} --fasta-input --threads $threads\
		--output ${out}/${sample}.krak $num

	fi
done

if [[ ! -s ${out}/${sample}.krak ]] ; then error "Something went wrong with running kraken... Exiting."; fi




echo " "
echo "########################################################################################################"
echo "########################          RUNNING KRAKEN-TRANSLATE ON OUTPUT            ########################"
echo "########################################################################################################"
echo " "

for file in ${out}/*.krak; do
	kraken-translate --db ${KRAKEN_DB} $file > ${file}en
	if [[ ! -s ${file}en ]] ; then error "Something went wrong with running kraken-translate... Exiting."; fi
	rm $file
done



echo " "
echo "########################################################################################################"
echo "########################            MAKING KRONAGRAM OF ALL FILES               ########################"
echo "########################################################################################################"
echo " "

#use custom script to summarize kraken file to krona format
for file in ${out}/*.kraken; do
	${SOFT}/kraken_to_krona.py $file > ${file%.*}.krona
	if [[ ! -s ${file%.*}.krona ]]; then error "Something went wrong with making krona file from kraken file. Exiting..."; fi
done

#use kronatools to make kronagrams of all samples in one html file
ktImportText -o ${out}/kronagram.html ${out}/*krona
if [[ ! -s ${out}/kronagram.html ]]; then error "Something went wrong with running KronaTools to make kronagram. Exiting..."; fi



echo " "
echo "########################################################################################################"
echo "########################         FINISHED RUNNING KRAKEN PIPELINE!!!            ########################"
echo "########################################################################################################"
echo " "






