#!/bin/bash -l

##############################################################################################################################################################
#
# This script is a modified pipeline from the program 'BLOBOLOGY', which produces a GC vs Coverage plot of the contigs, helping visualize bins and phylogenetic 
# composition.The original script has been modified to be better suited for running on clusters.
#
# Author of original pipeline: Sujai Kumar (https://github.com/blaxterlab). Author modifications: German Uritskiy. I do not take any credit for the original pipeline.
# For questions, bugs, and suggestions, contact German Uritskiy at guritsk1@jhu.edu.
# 
##############################################################################################################################################################



help_message () {
	echo "Usage: ./blobology.sh [options] -a assembly.fasta-1 reads_1.fastq -2 reads_2.fastq -o output_dir"
	echo "Options:"
	echo ""
	echo "	-a STR		assembly fasta file"
	echo "	-1 STR          forward fastq reads"
	echo "	-2 STR          reverse fastq reads" 
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads"
	echo "	-n INT		number of contigs to plot (default=ALL)"
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
threads=1; out="false"; read_1="false"; read_2="false"; n_contigs=1000000000; ASSEMBLY="false"

# Load in options
while getopts ht:a:o:1:2:n: option; do
	case "${option}" in
		h) help_message; exit 1;;
		t) threads=${OPTARG};;
		a) ASSEMBLY=$OPTARG;;
		o) out=${OPTARG};;
		1) reads_1=$OPTARG;;
		2) reads_2=$OPTARG;;
		n) n_contigs=$OPTARG;;
	esac
done



########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$reads_1" = "false" ] || [ "$reads_2" = "false" ] || [ "$threads" = "false" ] || [ "$ASSEMBLY" = "false" ]; then 
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure the meta-scripts folder path is correctly set in config.sh file"
fi

#  Checks for NCBI_nt database for BLAST
if [ ! -f "${BLASTDB}/nt.00.nhd" ]; then
	error "The file ${BLASTDB}/nt.00.nhd doesnt exist, which likely means that you havent set the correct path to your NCBI_nt database or\
	 havent downloaded it. Please look and the the script and manually set the correct path to NCBI_nt, and follow these steps to download the\
	  database: wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz; for a in nt.*.tar.gz; do tar xzf $a; done"
fi


#  Checks for NCBI_tax database for BLAST
if [ ! -f "${TAXDUMP}/citations.dmp" ]; then
	error "The file ${TAXDUMP}/citations.dmp doesnt exist, which likely means that you havent set the correct path to your NCBI_tax database,\
	 or havent downloaded it yet. Please look and the the script and manually set the correct path to NCBI_tax, and follow these steps to download \
	 the database: wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz; tar -xvf taxdump.tar.gz;"

fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################






echo " "
echo "########################################################################################################"
echo "########################       ASIGN TAXONOMY TO CONTIGS WITH MEGABLAST         ########################"
echo "########################################################################################################"
echo " "

mkdir $out
assembly=${ASSEMBLY##*/}
SAMPLE=${assembly%.*}

## the following command selects desired number of contigs at random to classify and plot:
comm "Choosing $n contigs from $ASSEMBLY"
${SOFT}/blobology/fastaqual_select.pl -f $ASSEMBLY -s r -n $n_contigs > ${out}/$assembly

comm "Running BLASTN"
blastn -task megablast -query ${out}/$assembly -db ${BLASTDB}/nt\
 -evalue 1e-5 -num_threads $threads -max_target_seqs 1 -outfmt '6 qseqid sseqid staxids'\
 | cut -f 1,3 > ${out}/${SAMPLE}.nt.1e-5.megablast

if [[ ! -s ${out}/${SAMPLE}.nt.1e-5.megablast ]] ; then 
	error "Something went wrong with asigning taxonomy to contigs with megablast. Exiting";
fi


echo " "
echo "########################################################################################################"
echo "########################          MAP READS TO ASSEMBLY WITH BOWTIE2            ########################"
echo "########################################################################################################"
echo " "

comm "Indexing ${out}/$assembly"
bowtie2-build -q ${out}/$assembly ${out}/$assembly

comm "Alligning $reads_1 and $$reads_2 to ${out}/$assembly with bowtie2"
${SOFT}/blobology/shuffleSequences_fastx.pl 4 <(cat $reads_1) <(cat $reads_2) > ${out}/tmp
bowtie2 -x ${out}/$assembly --very-fast-local -k 1 -t -p $threads --reorder --mm -U ${out}/tmp\
 | samtools view -S -b -T $threads > ${out}/${SAMPLE}.bowtie2.bam
rm ${out}/tmp

if [[ ! -s ${out}/${SAMPLE}.bowtie2.bam ]] ; then 
	error "Something went wrong with aligning reads back to the contigs with bowtie2. Exiting."; 
fi




echo " "
echo "########################################################################################################"
echo "########################     MAKE BLOB FILE FROM BLAST AND BOWTIE2 OUTPUTS      ########################"
echo "########################################################################################################"
echo " "


${SOFT}/blobology/gc_cov_annotate.pl --blasttaxid ${out}/${SAMPLE}.nt.1e-5.megablast\
 --assembly ${out}/$assembly --bam ${out}/${SAMPLE}.bowtie2.bam\
 --out ${out}/${SAMPLE}.blobplot --taxdump $TAXDUMP --taxlist species genus family order subclass phylum superkingdom
if [[ ! -s ${out}/${SAMPLE}.blobplot ]]; then 
	error "Something went wrong with making the blob file from the .bam and .megablast files. Exiting."; fi
comm "blobplot text file saved to ${out}/${SAMPLE}.blobplot"


echo " "
echo "########################################################################################################"
echo "########################              MAKE FINAL BLOBPLOT IMAGES                ########################"
echo "########################################################################################################"
echo " "

${SOFT}/blobology/makeblobplot.R ${out}/${SAMPLE}.blobplot 0.01 taxlevel_order
${SOFT}/blobology/makeblobplot.R ${out}/${SAMPLE}.blobplot 0.01 taxlevel_phylum
${SOFT}/blobology/makeblobplot.R ${out}/${SAMPLE}.blobplot 0.01 taxlevel_superkingdom
## The output file blobplot.txt can also be loaded into Blobsplorer - see github.com/mojones/blobsplorer
if [[ ! -s ${out}/${SAMPLE}.blobplot.taxlevel_phylum.png ]]; then 
	error "Something went wrong with making the plots from the blob file. Exiting."; 
fi


echo " "
echo "########################################################################################################"
echo "########################      BLOBPLOT PIPELINE FINISHED SUCCESSFULLY!!!        ########################"
echo "########################################################################################################"
echo " "


