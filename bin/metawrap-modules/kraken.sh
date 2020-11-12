#!/usr/bin/env bash

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
	echo ""
	echo "Run on any number of fasta assembly files and/or or paired-end reads."
	echo "Usage: metaWRAP kraken [options] -o output_dir assembly.fasta reads_1.fastq reads_2.fastq ..."
	echo "Options:"
	echo "" 
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads"
	echo "	-s INT		read subsampling number (default=all)"
	echo "	--no-preload	do not pre-load the kraken DB into memory (slower, but lower memory requirement)"
	echo ""
	echo "	Note: you may pass any number of sequence files with the following extensions:"
	echo "	*.fa *.fasta (assumed to be assembly files) or *_1.fastq and *_2.fastq (assumed to be paired)"
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

# Set defaults
threads=1; out="false"; depth="all"; preload=true


# load in params
OPTS=`getopt -o ht:o:s: --long help,no-preload -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
	case "$1" in
		-t) threads=$2; shift 2;;
		-o) out=$2; shift 2;;
		-s) depth=$2; shift 2;;
		-h | --help) help_message; exit 1; shift 1;;
		--no-preload) preload=false; shift 1;;
		--) help_message; exit 1; shift; break ;;
		*) break;;
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
	error "The folder $KRAKEN_DB doesnt exist. Please consult the metaWRAP database guite to download and build the KRAKEN database"
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################



########################################################################################################
########################              RUNNING KRAKEN ON ALL FILES               ########################
########################################################################################################
announcement "RUNNING KRAKEN ON ALL FILES"

# setting up the output folder
if [ ! -d $out ]; then 
	mkdir $out;
else 
	echo "Warning: $out already exists."
fi

# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	#process fastq files
	if [[ $num == *"_1.fastq" ]]; then
		reads_1=$num
		reads_2=${num%_*}_2.fastq
		if [ "$reads_2" = "false" ]; then error "$reads_2 does not exist. Exiting..."; fi
		
		tmp=${reads_1##*/}
		sample=${tmp%_*}
		comm "Now processing $reads_1 and $reads_2 with $threads threads"

		# if sampling depth is specified, randomly subsample the fastq reads
		if [ ! "$depth" = "all" ]; then
			comm "subsampling down to $depth reads..." 
			paste $reads_1 $reads_2 | \
			 awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | \
			 shuf | head -n $depth | sed 's/\t\t/\n/g' | \
			 awk -F"\t" '{print $1 > "'"${out}/tmp_1.fastq"'"; print $2 > "'"${out}/tmp_2.fastq"'"}'
			reads_1=${out}/tmp_1.fastq
			reads_2=${out}/tmp_2.fastq
			comm "Subsampling done. Starting KRAKEN..."
			if [ ! -s $reads_1 ]; then error "something went wrong with subsampling sequences. Exiting..."; fi
		fi

		if [ ! -s $reads_1 ]; then error "$reads_1 doesnt exist. Exiting..."; fi

		if [ "$preload" = true ]; then
			kraken --db ${KRAKEN_DB} --fastq-input --paired --threads $threads --preload --output ${out}/${sample}.krak $reads_1 $reads_2
		else
			kraken --db ${KRAKEN_DB} --fastq-input --paired --threads $threads --output ${out}/${sample}.krak $reads_1 $reads_2
		fi
	
		if [[ $? -ne 0 ]] || [[ ! -s ${out}/${sample}.krak ]]; then error "Something went wrong with running kraken on $reads_1 and $reads_2 . Exiting..."; fi
			
		if [ ! "$depth" = "all" ]; then rm ${out}/tmp_1.fastq ${out}/tmp_2.fastq; fi
	fi

	#process fasta files
	if [[ $num == *"fa" ]] || [[ $num == *"fasta" ]]; then
		tmp=${num##*/}
		sample=${tmp%.*}
		comm "Now processing $num with $threads threads"
		
		if [ "$preload" = true ]; then
			kraken --db ${KRAKEN_DB} --fasta-input --threads $threads --preload --output ${out}/${sample}.krak $num
		else
			kraken --db ${KRAKEN_DB} --fasta-input --threads $threads --output ${out}/${sample}.krak $num
		fi		

		if [[ $? -ne 0 ]] || [[ ! -s ${out}/${sample}.krak ]]; then error "Something went wrong with running kraken on ${num}. Exiting..."; fi
	fi
done

# check if any files were processed
if [[ $( ls $out | grep ".krak" | wc -l ) -eq 0 ]]; then 
	comm "No fasta or fastq files detected! (must be in .fastq .fa .fastq or .fq format)"
	help_message; exit 1
fi


########################################################################################################
########################          RUNNING KRAKEN-TRANSLATE ON OUTPUT            ########################
########################################################################################################
announcement "RUNNING KRAKEN-TRANSLATE ON OUTPUT"

for file in ${out}/*.krak; do
	comm "Translating $file"
	kraken-translate --db ${KRAKEN_DB} $file > ${file}en
	if [[ ! -s ${file}en ]] ; then error "Something went wrong with running kraken-translate... Exiting."; fi
	rm $file
done




########################################################################################################
########################            MAKING KRONAGRAM OF ALL FILES               ########################
########################################################################################################
announcement "MAKING KRONAGRAM OF ALL FILES"

#use custom script to summarize kraken file to krona format
for file in ${out}/*.kraken; do
	${SOFT}/kraken_to_krona.py $file > ${file%.*}.krona
	if [[ ! -s ${file%.*}.krona ]]; then error "Something went wrong with making krona file from kraken file. Exiting..."; fi
done

#use kronatools to make kronagrams of all samples in one html file
ktImportText -o ${out}/kronagram.html ${out}/*krona
if [[ ! -s ${out}/kronagram.html ]]; then error "Something went wrong with running KronaTools to make kronagram. Exiting..."; fi


########################################################################################################
########################         FINISHED RUNNING KRAKEN PIPELINE!!!            ########################
########################################################################################################
announcement "FINISHED RUNNING KRAKEN PIPELINE!!!"

