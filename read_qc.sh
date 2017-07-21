#!/usr/bin/env bash

###########################################################################################################################################################
#       	                                                                                                                                          #
# This script is meant to be a comprehensive solution to QC new HiSeq reads in preparation for assembly and other operations.                             #
# The main things this pipeline accomplishes are read trimming based on quality scores, and removal of human sequences.                                   #
# The script also produces a FASTQC report before and after the procedures.                                                                               #
#	                                                                                                                                 		  #
###########################################################################################################################################################

help_message () {
	echo "Usage: ./read_qc.sh [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir"
	echo "Options:"
	echo ""
	echo "	-1 STR          forward fastq reads"
	echo "	-2 STR          reverse fastq reads" 
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads"
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

thread=sfalse; mem=false; out="false"; reads_1="false"; reads_2="false"
# Load in options
while getopts ht:o:1:2: option; do
        case "${option}" in
		h) help_message; exit 1;;
		t) threads=${OPTARG};;
		o) out=${OPTARG};;
		1) reads_1=$OPTARG;;
		2) reads_2=$OPTARG;;
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

if [ ! -s ${BMTAGGER_DB}/hg38.bitmask ]; then
	error "${BMTAGGER_DB}/hg38.bitmask file doesnt exist. Please configure your bmtagger human genome index"
fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################




mkdir $out

echo " "
echo "########################################################################################################"
echo "########################                 MAKING PRE-QC REPORT                   ########################"
echo "########################################################################################################"
echo " "

mkdir ${out}/pre-QC_report
fastqc -q -t $threads -o ${out}/pre-QC_report -f fastq $reads_1 $reads_2
rm ${out}/pre-QC_report/*zip

if [[ ! -s ${out}/pre-QC_report/${reads_1%.*}_fastqc.html ]]; then error "Something went wrong with making pre-QC fastqc report. Exiting."; fi
comm "pre-qc report saved to: ${out}/pre-QC_report"


echo " "
echo "########################################################################################################"
echo "########################                 RUNNING TRIM-GALORE                    ########################"
echo "########################################################################################################"
echo " "

trim_galore --no_report_file --paired -o $out $reads_1 $reads_2

tmp=${reads_1%_*}; sample=${tmp##*/}

# Fix the naming of the trimmed reads files:
mv ${out}/${sample}_1_val_1.fq ${out}/${sample}_1.trimmed.fastq
mv ${out}/${sample}_2_val_2.fq ${out}/${sample}_2.trimmed.fastq
if [[ ! -s ${out}/${sample}_1.trimmed.fastq ]]; then error "Something went wrong with trimming the reads. Exiting."; fi
comm "Trimmed reads saved to: ${out}/${sample}*.trimmed.fastq"


echo " "
echo "########################################################################################################"
echo "########################               REMOVING HUMAN SEQUENCES                 ########################"
echo "########################################################################################################"
echo " "

mkdir ${out}/bmtagger_tmp
bmtagger.sh -b ${BMTAGGER_DB}/hg38.bitmask -x ${BMTAGGER_DB}/hg38.srprism -T ${out}/bmtagger_tmp -q1\
 -1 ${out}/${sample}_1.trimmed.fastq\
 -2 ${out}/${sample}_2.trimmed.fastq\
 -o ${out}/${sample}.bmtagger.list
if [[ ! -s ${out}/${sample}.bmtagger.list ]]; then warning "Something went wrong with finding human reads (unlikely that there are none...). Just sayin'"; fi


comm "Now sorting out found human reads from the main fastq files..."
${SOFT}/skip_human_reads.py ${out}/${sample}.bmtagger.list ${out}/${sample}_1.trimmed.fastq > ${out}/${sample}_1.clean.fastq
${SOFT}/skip_human_reads.py ${out}/${sample}.bmtagger.list ${out}/${sample}_2.trimmed.fastq > ${out}/${sample}_2.clean.fastq

comm "Now sorting out found human reads and putting them into a new file... for science..."
${SOFT}/select_human_reads.py ${out}/${sample}.bmtagger.list ${out}/${sample}_1.trimmed.fastq > ${out}/human_reads_1.fastq
${SOFT}/select_human_reads.py ${out}/${sample}.bmtagger.list ${out}/${sample}_2.trimmed.fastq > ${out}/human_reads_2.fastq
if [[ ! -s ${out}/${sample}_1.clean.fastq ]]; then error "Something went wrong with removing human reads with bmtagger. Exiting."; fi


mv ${out}/${sample}_1.clean.fastq ${out}/final_pure_reads_1.fastq
mv ${out}/${sample}_2.clean.fastq ${out}/final_pure_reads_2.fastq
comm "Contamination-free and trimmed reads are stored in: ${out}/final_pure_reads_1.fastq and ${out}/final_pure_reads_2.fastq"



echo " "
echo "########################################################################################################"
echo "########################                 MAKING POST-QC REPORT                  ########################"
echo "########################################################################################################"
echo " "

mkdir ${out}/post-QC_report
fastqc -t $threads -o ${out}/post-QC_report -f fastq ${out}/final_pure_reads_1.fastq and ${out}/final_pure_reads_2.fastq
rm ${out}/post-QC_report/*zip

if [[ ! -s ${out}/post-QC_report/final_pure_reads_1_fastqc.html ]]; then error "Something went wrong with making post-QC fastqc report. Exiting."; fi
comm "pre-qc report saved to: ${out}/post-QC_report"


# Remove intermediate files
rm -r ${out}/bmtagger_tmp


echo " "
echo "########################################################################################################"
echo "########################              READ QC PIPELINE COMPLETE!!!              ########################"
echo "########################################################################################################"
echo " "








