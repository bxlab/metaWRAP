#!/bin/bash -l

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
	echo "Usage: ./assembly.sh [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir -t THREADS -m MEMORY"
	echo "Options:"
	echo ""
	echo "	-1 STR          forward fastq reads"
	echo "	-2 STR          reverse fastq reads" 
	echo "	-o STR          output directory"
	echo "	-m INT          memory in GB"
	echo "	-t INT          number of threads"
	echo "";}
# function to print out error messages
error () {
	echo ""; echo "************************************************************************************************************************"
	echo "*****************************************            ERROR!           **************************************************"
	echo $1
	echo "************************************************************************************************************************"; echo ""; exit 1; }
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


threads=1; out="false"; read_1="false"; read_2="false"
# Load in options
while getopts ht:m:o:1:2: option; do
	case "${option}" in
		h) help_message; exit 1;;
		t) threads=${OPTARG};;
		m) mem=$OPTARG;;
		o) out=${OPTARG};;
		1) reads_1=$OPTARG;;
		2) reads_2=$OPTARG;;
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





echo " "
echo "########################################################################################################"
echo "########################               ASSEMBLING WITH METASPADES               ########################"
echo "########################################################################################################"
echo " "

mkdir $out

comm "Using reads $reads_1 and $reads_2 for assembly,"
if [ ! -f "$reads_1" ] || [ ! -f "$reads_2" ]; then error "Read files $reads_1 and/or $reads_2 dont exist. Exiting."; fi

metaspades.py --tmp-dir ${out}/metaspades.tmp -t $threads -m $mem -o ${out}/metaspades -1 $reads_1 -2 $reads_2 
rm -r ${out}/metaspades.tmp
if [ ! -f "${out}/metaspades/scaffolds.fasta" ]; then error "Something went wrong with metaSPAdes assembly. Exiting."; fi



echo " "
echo "########################################################################################################"
echo "########################              SORTING OUT UNASSEMBLED READS             ########################"
echo "########################################################################################################"
echo " "

#take only the scaffolds over 1.5kb (the minimum needed for binning with metaBAT):
${SOFT}/rm_short_contigs.py 1500 ${out}/metaspades/scaffolds.fasta > ${out}/metaspades/long_scaffolds.fasta
if [ ! -s ${out}/metaspades/long_scaffolds.fasta ]; then
	error "metaSPAdes produced no contigs over 1.5kb, so the rest of the pipeline wont work. Exiting."; fi


bwa index ${out}/metaspades/long_scaffolds.fasta

#sort out and store reads that dont map back to the assembly:
bwa mem -t $threads ${out}/metaspades/long_scaffolds.fasta $reads_1 $reads_2 | grep -v NM:i: > ${out}/unused_by_metaspades.sam
${SOFT}/sam_to_fastq.py ${out}/unused_by_metaspades.sam > ${out}/unused_by_metaspades.fastq
#rm ${out}/unused_by_metaspades.sam

if [[ ! -s ${out}/unused_by_metaspades.fastq ]]; then error "Something went wrong with pulling out unassembled reads. Exiting."; fi



echo " "
echo "########################################################################################################"
echo "########################        ASSEMBLING REMAINDER READS WITH MEGAHIT         ########################"
echo "########################################################################################################"
echo " "

rm -r ${out}/megahit
megahit -r ${out}/unused_by_metaspades.fastq -o ${out}/megahit -t $threads -m ${mem}000000000
if [ ! -f "${out}/megahit/final.contigs.fa" ]; then error "Something went wrong with reassembling with Megahit. Exiting."; fi
mv ${out}/unused_by_metaspades.fastq ${out}/metaspades/


echo " "
echo "########################################################################################################"
echo "########################         COMBINE AND FORMAT THE TWO ASSEMBLIES          ########################"
echo "########################################################################################################"
echo " "

${SOFT}/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa > ${out}/megahit/fixed.contigs.fa
${SOFT}/rm_short_contigs.py 500 ${out}/megahit/fixed.contigs.fa > ${out}/megahit/long.contigs.fa

cp ${out}/metaspades/long_scaffolds.fasta ${out}/combined_assembly.fasta
cat ${out}/megahit/long.contigs.fa >> ${out}/combined_assembly.fasta
${SOFT}/sort_contigs.py ${out}/combined_assembly.fasta > ${out}/final_assembly.fasta
if [[ ! -s ${out}/final_assembly.fasta ]]; then error "Something went wrong with joining the assemblies. Exiting."; fi
rm ${out}/combined_assembly.fasta


echo " "
echo "########################################################################################################"
echo "########################             RUNNING ASSEMBLY QC WITH QUAST             ########################"
echo "########################################################################################################"
echo " "

cp ${out}/metaspades/scaffolds.fasta ${out}/metaspades/metaspades_assembly.fasta
quast -t $threads -o ${out}/QUAST_out -m 500 ${out}/final_assembly.fasta ${out}/metaspades/metaspades_assembly.fasta
cp ${out}/QUAST_out/report.html ${out}/assembly_report.html

if [[ ! -s ${out}/assembly_report.html ]]; then error "Something went wrong with running QUAST. Exiting."; fi


echo " "
echo "########################################################################################################"
echo "########################      ASSEMBLY PIPELINE COMPLETED SUCCESSFULLY!!!       ########################"
echo "########################################################################################################"
echo " "





