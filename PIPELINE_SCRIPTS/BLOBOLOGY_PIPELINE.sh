#!/bin/bash -l
#SBATCH
#SBATCH --job-name=BLOBOLOGY
#SBATCH --partition=lrgmem
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=900G

echo " "
echo "##############################################################################################################################################################"
echo "#"
echo "# INPUT: sbatch BLOBOLOGY_PIPELINE.sh assembly.fa reads_1.fastq reads_2.fastq"
echo "# This script is a modified pipeline from the program 'BLOBOLOGY'. The original script has been modified to be better suited for running on clusters."
echo "#"
echo "# Author of original pipeline: Sujai Kumar (https://github.com/blaxterlab). Author modifications: German Uritskiy. I do not take any credit for the original pipeline."
echo "# For questions, bugs, and suggestions, contact German Uritskiy at guritsk1@jhu.edu."
echo "# "
echo "##############################################################################################################################################################"
echo " "


####==========================================================================
#### GENERAL CONFIG VARIABLES
####==========================================================================

## Location of local installation of nt blast database
## (not needed if using blast remotely, which is slower).
## The NCBI nt databases can be downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/blast/db/ using the following command:
# wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
# for a in nt.*.tar.gz; do tar xzf $a; done
BLASTDB=/home-2/guritsk1@jhu.edu/scratch/RefSeq/NCBI_nt

## Location of NCBI tar gunzipped directory downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
## Default: current directory
TAXDUMP=/home-2/guritsk1@jhu.edu/scratch/RefSeq/NCBI_tax


# Don't forget to set the correct path to the SOFTWARE folder, containg all the needed packages. This must be an absolute path:"
SOFT=/home-2/guritsk1@jhu.edu/scratch/METAGENOME_PIPELINE/SOFTWARE/

threads=48
ASSEMBLY=$1
SAMPLE=${ASSEMBLY##*/}
reads_1=$2
reads_2=$3
out=${SAMPLE%.*}.blobology

# Checks for correct script input
if [ "$#" -lt 3 ]; then
echo "************************************************************************************************************************"
echo "*****************************************    ERROR: INVALID INPUT     **************************************************"
echo " "
echo "Please run this script like this: ./BLOBOLOGY_PIPELINE assembly.fa sample_A.fastq sampleA_2.fastq"
echo " "
echo "************************************************************************************************************************"
exit 1
fi

# Checks for correctly configures SOFTWARE folder
if [ ! -d "$SOFT" ]; then
echo "************************************************************************************************************************"
echo "*******************************    ERROR: INVALID PATH TO SOFTWARE FOLDER:   *******************************************"
echo " "
echo "The folder $SOFT doesnt exist. Please look and the the script and manually set the correct path to the 'SOFTWARE' folder"
echo " "
echo "************************************************************************************************************************"
exit 1
fi


#  Checks for NCBI_nt database for BLAST
if [ ! -f "${BLASTDB}/nt.00.nhd" ]; then
echo "************************************************************************************************************************"
echo "*******************************      ERROR: NCBI_NT DATABASE NOT FOUND       *******************************************"
echo " "
echo "The file ${BLASTDB}/nt.00.nhd  doesnt exist, which likely means that you havent set the correct path to your NCBI_nt database. "
echo "Please look and the the script and manually set the correct path to NCBI_nt, or follow these steps to download the database:"
echo " "
echo "The NCBI nt databases can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/db/ using the following command:"
echo 'wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"'
echo "for a in nt.*.tar.gz; do tar xzf $a; done"
echo " "
echo "************************************************************************************************************************"
exit 1
fi


#  Checks for NCBI_tax database for BLAST
if [ ! -f "${TAXDUMP}/citations.dmp" ]; then
echo "************************************************************************************************************************"
echo "*******************************      ERROR: NCBI_TAX DATABASE NOT FOUND      *******************************************"
echo " "
echo "The file ${TAXDUMP}/citations.dmp  doesnt exist, which likely means that you havent set the correct path to your NCBI_tax database. "
echo "Please look and the the script and manually set the correct path to NCBI_tax, or follow these steps to download the database:"
echo " "
echo "The NCBI tax databases can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/db/ using the following command:"
echo "wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
echo "tar -xvf taxdump.tar.gz"
echo " "
echo "************************************************************************************************************************"
exit 1
fi





echo " "
echo "########################################################################################################"
echo "########################       ASIGN TAXONOMY TO CONTIGS WITH MEGABLAST         ########################"
echo "########################################################################################################"
echo " "

mkdir $out
## the following command selects 10000 contigs at random:
echo "Choosing contigs"
RND_ASSEMBLY=random10k.$SAMPLE
#cp $ASSEMBLY ${out}/$RND_ASSEMBLY
${SOFT}/blobology/fastaqual_select.pl -f $ASSEMBLY -s r -n 10000 > ${out}/$RND_ASSEMBLY

echo "Running BLASTN"
blastn -task megablast -query ${out}/$RND_ASSEMBLY -db ${BLASTDB}/nt\
 -evalue 1e-5 -num_threads $threads -max_target_seqs 1 -outfmt '6 qseqid sseqid staxids'\
 | cut -f 1,3 > ${out}/random10k.assembly.nt.1e-5.megablast

if [[ ! -s ${out}/random10k.assembly.nt.1e-5.megablast ]] ; then echo "Something went wrong with asigning taxonomy to contigs with megablast. Exiting"; exit 1; fi


echo " "
echo "########################################################################################################"
echo "########################          MAP READS TO ASSEMBLY WITH BOWTIE2            ########################"
echo "########################################################################################################"
echo " "

echo "Indexing ${out}/$RND_ASSEMBLY"
bowtie2-build ${out}/$RND_ASSEMBLY ${out}/$RND_ASSEMBLY

echo "Alligning $reads_1 and $$reads_2 to ${out}/$RND_ASSEMBLY with bowtie2"
${SOFT}/blobology/shuffleSequences_fastx.pl 4 <(cat $reads_1) <(cat $reads_2) > ${out}/tmp
bowtie2 -x ${out}/$RND_ASSEMBLY --very-fast-local -k 1 -t -p $threads --reorder --mm -U ${out}/tmp\
 | samtools view -S -b -T $threads > ${out}/${RND_ASSEMBLY%.*}.bowtie2.bam
rm ${out}/tmp

if [[ ! -s ${out}/${RND_ASSEMBLY%.*}.bowtie2.bam ]] ; then echo "Something went wrong with alignning reads back to the contigs with bowtie2. Exiting."; exit 1; fi

# Note: shuffleSequences_fastx.pl works on gunzipped files, hence the need to zcat inline



echo " "
echo "########################################################################################################"
echo "########################     MAKE BLOB FILE FROM BLAST AND BOWTIE2 OUTPUTS      ########################"
echo "########################################################################################################"
echo " "

${SOFT}/blobology/gc_cov_annotate.pl --blasttaxid ${out}/random10k.assembly.nt.1e-5.megablast\
 --assembly ${out}/$RND_ASSEMBLY --bam ${out}/${RND_ASSEMBLY%.*}.bowtie2.bam\
 --out ${out}/${ASSEMBLY%.*}.blobplot --taxdump $TAXDUMP --taxlist species order phylum superkingdom
if [[ ! -s ${out}/${ASSEMBLY%.*}.blobplot ]]; then echo "Something went wrong with making the blob file from the .bam and .megablast files. Exiting."; exit 1; fi


echo " "
echo "########################################################################################################"
echo "########################              MAKE FINAL BLOBPLOT IMAGES                ########################"
echo "########################################################################################################"
echo " "

${SOFT}/blobology/makeblobplot.R ${out}/${ASSEMBLY%.*}.blobplot 0.01 taxlevel_order
${SOFT}/blobology/makeblobplot.R ${out}/${ASSEMBLY%.*}.blobplot 0.01 taxlevel_phylum
${SOFT}/blobology/makeblobplot.R ${out}/${ASSEMBLY%.*}.blobplot 0.01 taxlevel_superkingdom
## The output file blobplot.txt can also be loaded into Blobsplorer - see github.com/mojones/blobsplorer
if [[ ! -s ${out}/${ASSEMBLY%.*}.blobplot.taxlevel_phylum.png ]]; then echo "Something went wrong with making the plots from the blob file. Exiting."; exit 1; fi


echo " "
echo "########################################################################################################"
echo "########################             BLOBPLOT PIPELINE FINISHED!!!              ########################"
echo "########################################################################################################"
echo " "


