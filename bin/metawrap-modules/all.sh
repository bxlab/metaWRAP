#!/usr/bin/env bash

##############################################################################################################################################################
# Master metaWRAP script that calls the read_qc, assembly, blobology, binning, and kraken modules
##############################################################################################################################################################
config_file=$(which config-metawrap)
source $config_file

help_message () {
	echo "";
	echo "Usage: ./metaWRAP all [options] -o output_folder -1 raw_readsA_1.fastq -2 raw_readsA_2.fastq"
	echo "Note: currently only works on one pair of reads. To run run multiple samples. However individual modules can handle more samples."
	echo "Options:"
	echo ""
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT          memory in GB (default=4)"
	echo "	-1 STR		forward read file"
	echo "	-2 STR		reverse read file"
	echo "	--fast		skips some non-essential and computationally expensive parts of the pipeline:"
	echo "				-skips human sequence removal and pre-qc fastq report in read_qc module"
	echo "  	                        -replaces metaspades with megahit as primary assembler"
	echo "          	                -subsamples Kraken module to 10k reads instead of 1M reads"
	echo "                  	        -skips bin reassembly in Binning module, and skips Kraken on bins"
	echo "                          	-subsammples to 10k contigs and skips bin annotation in Blobology module"
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "!"; }


########################################################################################################
########################                  LOADING IN PARAMETERS                 ########################
########################################################################################################

announcement "Running 'ALL' modules!"

#defaults:
threads=1; mem=4; out="false"; reads_1="false"; reads_2="false"; fast=false

OPTS=`getopt -o ht:m:o:1:2: --long help,fast -- "$@"`
# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
		-m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
                -1) reads_1=$2; shift 2;;
                -2) reads_2=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
                --fast) fast=true; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################
comm "checking if all parameters are entered"

if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ] ; then
	comm "Error when checking -1 -2 or -o options"
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure the scripts folder is in the same folder as the other reads_qc"
fi

if [ -d $out ]; then warning "WARNING: Directory $out already exists..."; fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################

mkdir $out
################################### READ_QC MODULE #######################################
announcement "STARTING READ_QC MODULE"
if [ "$fast" = true ]; then 
	comm "Skipping pre-trimming FASTQC report and human sequence removing with bmtagger (because --fast option was specified)"
	cmd="${PIPES}/read_qc.sh -t $threads -1 $reads_1 -2 $reads_2 -o ${out}/read_qc_out --skip-bmtagger --skip-pre-qc-report"
else 
	cmd="${PIPES}/read_qc.sh -t $threads -1 $reads_1 -2 $reads_2 -o ${out}/read_qc_out"
fi
echo $cmd
$cmd

if [ ! -s ${out}/read_qc_out/final_pure_reads_1.fastq ]; then
	error "read_qc module did not complete successfully. Aborting pipeline.";
else
	comm "read_qc module finished successfully!"
fi

################################### ASSEMBLY MODULE #######################################
announcement "STARTING ASSEMBLY MODULE"
if [ "$fast" = true ]; then
	comm "Skipping metaspades assembly, going straight to megahit assembly (because --fast option was specified)"
	cmd="${PIPES}/assembly.sh -t $threads -m $mem\
 	-1 ${out}/read_qc_out/final_pure_reads_1.fastq\
 	-2 ${out}/read_qc_out/final_pure_reads_2.fastq\
 	-o $out/assembly_out"
else
	cmd="${PIPES}/assembly.sh -t $threads -m $mem\
	-1 ${out}/read_qc_out/final_pure_reads_1.fastq\
	-2 ${out}/read_qc_out/final_pure_reads_2.fastq\
	-o $out/assembly_out --metaspades"
fi
echo $cmd
$cmd

if [ ! -s ${out}/assembly_out/final_assembly.fasta ]; then
	error "Assembly module did not complete successfully. Aborting pipeline.";
else
	comm "Assembly module completed successfully!"
fi



#################################### KRAKEN MODULE #######################################
announcement "STARTING KRAKEN MODULE"
if [ "$fast" = true ]; then
	comm "Subsampling to 10k reads and ignoring assembly (because --fast option was specified)"
	cmd="${PIPES}/kraken.sh -t $threads -o ${out}/kraken_out -s 10000\
	 ${out}/read_qc_out/final_pure_reads_1.fastq\
	 ${out}/read_qc_out/final_pure_reads_2.fastq"
else
        cmd="${PIPES}/kraken.sh -t $threads -o ${out}/kraken_out -s 1000000\
         ${out}/assembly_out/final_assembly.fasta\
         ${out}/read_qc_out/final_pure_reads_1.fastq\
         ${out}/read_qc_out/final_pure_reads_2.fastq"
fi
echo $cmd
$cmd

if [ -s ${out}/kraken_out/final_assembly.krona ] || [ -s ${out}/kraken_out/final_pure_reads.krona ] ; then
	comm "KRAKEN module completed sucessfully!"
else
	error "KRAKEN module did not complete successfully. Aborting pipeline.";
fi


################################### BINNING MODULE #######################################
announcement "STARTING BINNING MODULE"
if [ "$fast" = true ]; then
	comm "Skipping reassembly and selecting good bins (because --fast option was specified)"
	cmd="${PIPES}/binning.sh -t $threads -m $mem --skip-reassembly\
	 -a $out/assembly_out/final_assembly.fasta -o ${out}/binning_out\
	 $out/read_qc_out/final_pure_reads_1.fastq $out/read_qc_out/final_pure_reads_2.fastq"
else
	cmd="${PIPES}/binning.sh -t $threads -m $mem --checkm-good-bins --checkm-best-bins\
         -a $out/assembly_out/final_assembly.fasta -o ${out}/binning_out\
         $out/read_qc_out/final_pure_reads_1.fastq $out/read_qc_out/final_pure_reads_2.fastq"
fi
echo $cmd
$cmd

if [ ! -s ${out}/binning_out/metabat2_bins/bin.unbinned.fa ]; then
	error "Binning module did not complete successfully. Aborting pipeline."
else
	comm "Binning module completed successfully!"
fi


################################### BLOBOLOGY MODULE ######################################
announcement "STARTING BLOBOLOGY MODULE"
if [ "$fast" = true ]; then
	comm "Subsampling to 10k contigs and not annotating bins (because --fast option was speciefied)"
	cmd="${PIPES}/blobology.sh -t $threads -a $out/assembly_out/final_assembly.fasta\
	 --subsample 10000\ 
	 -1 $out/read_qc_out/final_pure_reads_1.fastq\
	 -2 $out/read_qc_out/final_pure_reads_2.fastq\
	 -o ${out}/blobology_out"
else
	cmd="${PIPES}/blobology.sh -t $threads -a $out/assembly_out/final_assembly.fasta\
         --bins ${out}/binning_out/good_bins\
         -1 $out/read_qc_out/final_pure_reads_1.fastq\
         -2 $out/read_qc_out/final_pure_reads_2.fastq\
         -o ${out}/blobology_out"
fi
echo $cmd
$cmd

if [ ! -s ${out}/blobology_out/final_assembly.blobplot ]; then
	error "Blobology module did not complete successfully. Aborting pipeline."
else
	comm "Blobology module completed successfully!"
fi

################################ KRAKEN MODULE ON BINS  ####################################
if [ ! "$fast" = true ]; then
	announcement "STARTING KRAKEN MODULE ON BINS"
	cmd="${PIPES}/kraken.sh -t $threads -o ${out}/binning_out/kraken\
	 ${out}/binning_out/metabat2_bins/*.fa"
	echo $cmd
	$cmd
	
	if [ ! -s ${out}/binning_out/kraken/bin.1.krona ]; then
		error "Kraken module on bins did not complete successfully. Aborting pipeline.";
	else
		comm "Kraken module completed sucessfully!"
	fi
else
	comm "Skipping runninng KRAKEN on bins (because --fast options specified)"
fi

######################## COPYING OVER ALL FIGURES AND REPORTS  #############################
announcement "Copying all figures and reports into 'metaWRAP_figures' folder for easy viewing"
mkdir ${out}/metaWRAP_figures
cp ${out}/read_qc_out/post-QC_report/final_pure_reads_1_fastqc.html ${out}/metaWRAP_figures/reads_1.fastqc.html
cp ${out}/read_qc_out/post-QC_report/final_pure_reads_2_fastqc.html ${out}/metaWRAP_figures/reads_2.fastqc.html
cp ${out}/assembly_out/assembly_report.html ${out}/metaWRAP_figures/
cp ${out}/binning_out/*.png ${out}/metaWRAP_figures/
cp ${out}/binning_out/*.stats ${out}/metaWRAP_figures/
cp ${out}/blobology_out/*.png ${out}/metaWRAP_figures/
cp ${out}//binning_out/insert_sizes.txt ${out}/metaWRAP_figures/

if [ ! "$fast" = true ]; then
	cp ${out}/read_qc_out/pre-QC_report/*_1_fastqc.html ${out}/metaWRAP_figures/raw_reads_1.fastqc.html
	cp ${out}/read_qc_out/pre-QC_report/*_2_fastqc.html ${out}/metaWRAP_figures/raw_reads_2.fastqc.html
fi

################################### END OF PIPELINE  #######################################
announcement "METAWRAP FINISHED SUCCESSFULLY! END OF PIPELINE."


