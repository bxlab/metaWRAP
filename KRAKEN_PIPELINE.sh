#!/bin/bash -l
#SBATCH
#SBATCH --job-name=KRAKEN_pipeline
#SBATCH --partition=shared
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G

echo " "
echo "##############################################################################################################################################################"
echo "#"
echo "# READS INPUT: sbatch KRAKEN_PIPELINE.sh output_folder sampleA_1.fastq sampleA_2.fastq [ sampleX_1.fastq sampleX_2.fastq ]"
echo "# ASSEMBLY INPUT: sbatch KRAKEN_PIPELINE.sh output_flder sampleA.fasta sampleA.fa"
echo "#"
echo "# This script is meant to be run on paired end reads (with extensions *_1.fastq and *_2.fastq) or assembled contigs (*.fa or *.fasta)."
echo "# The script runs KRAKEN on the sequences, then translates them to taxonomy form with kraken-translate. Then in-house scripts are used to "
echo "# parse out the taxonomy into the format for KRONA-TOOLS, colapse the file to save memory, and finally produce a prety kronagram with all the files."
echo "#"
echo "# NOTE: KRAKEN KronaTools requires instalation, and be sure the configure the right path to the KRAKEN folder and the KRAKEN database."
echo "#"
echo "# Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses."
echo "# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu."
echo "# "
echo "##############################################################################################################################################################"
echo " "




# Don't forget to set the correct path to the SOFTWARE folder, containg all the needed packages. This must be an absolute path:"
SOFT="/home-2/guritsk1@jhu.edu/scratch/METAGENOME_PIPELINE/SOFTWARE"
KRAKEN="/home-2/guritsk1@jhu.edu/scratch/Programs/kraken"
KRAKEN_DB="/home-2/guritsk1@jhu.edu/scratch/Programs/kraken/kraken_standard_database"
threads=24
mem=120



# Checks for correct script input
if [ "$#" -lt 2 ] || [[ $1 == *".fastq"* ]] || [[ $1 == *".fq"* ]] || [[ $1 == *".fasta"* ]] || [[ $1 == *".fa"* ]]; then
echo "************************************************************************************************************************"
echo "*****************************************    ERROR: INVALID INPUT     **************************************************"
echo " "
echo "Please enter input files in the following format:"
echo "READS INPUT: sbatch KRAKEN_PIPELINE.sh output_folder sampleA_1.fastq sampleA_2.fastq [ sampleX_1.fastq sampleX_2.fastq ]"
echo "or"
echo "ASSEMBLY INPUT: sbatch KRAKEN_PIPELINE.sh output_folder sampleA.fasta sampleX.fa"
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

#  Checks for KRAKEN instalation
if [ ! -f "${KRAKEN}/kraken" ]; then
echo "************************************************************************************************************************"
echo "*******************************       ERROR: KRAKEN INSTALATION NOT FOUND       *******************************************"
echo " "
echo "The file ${KRAKEN}/kraken doesnt exist, which means that you havent set the correct path to your KRAKEN instalation. "
echo "Please look and the the script and manually set the correct path to KRAKEN"
echo " "
echo "************************************************************************************************************************"
exit 1
fi

#Checks for present KREKEN database 
if [ ! -d "$KRAKEN_DB" ]; then
echo "************************************************************************************************************************"
echo "*******************************       ERROR: KRAKEN DATABASE NOT FOUND       *******************************************"
echo " "
echo "The folder $KRAKEN_DB doesnt exist. Please look and the the script and manually set the correct path to the KRAKEN standard database"
echo "If you do not have it yet, you can download it like this: "
echo "${KRAKEN}/kraken-build --standard --threads $threads --db ${KRAKEN}/kraken/KRAKEN_DB"
echo " "
echo "************************************************************************************************************************"
exit 1
fi










echo " "
echo "########################################################################################################"
echo "########################              RUNNING KRAKEN ON ALL FILES               ########################"
echo "########################################################################################################"
echo " "

# setting up the output folder
out=$1
mkdir $out

# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	if [[ $num == *"_1."*"fastq"* ]]; then 
		reads_1=$num
		reads_2="${reads_1/_1/_2}"
		tmp=${reads_1##*/}
		sample=${tmp%_*}
		echo "####       Now processing $reads_1 and $reads_2         ####"

		echo "${KRAKEN}/kraken --db ${KRAKEN_DB} --fastq-input --paired --threads $threads\
                 --output ${out}/${sample}.krak $reads_1 $reads_2"
	
		${KRAKEN}/kraken --db ${KRAKEN_DB} --fastq-input --paired --threads $threads\
		 --output ${out}/${sample}.krak $reads_1 $reads_2
	fi

	if [[ $num == *"fa"* ]] || [[ $num == *"fasta"* ]]; then
                reads=$num
                tmp=${reads##*/}
                sample=${tmp%_*}
                echo "####       Now processing $reads        ####"
		
		echo "${KRAKEN}/kraken --db ${KRAKEN_DB} --fasta-inputt --threads $threads\
                 --output ${out}/${sample}.krak $reads"
	
                ${KRAKEN}/kraken --db ${KRAKEN_DB} --fasta-input --threads $threads\
                 --output ${out}/${sample}.krak $reads
        fi
done

if [[ ! -s ${out}/${sample}.krak ]] ; then echo "Something went wrong with running kraken... Exiting."; exit 1; fi




echo " "
echo "########################################################################################################"
echo "########################          RUNNING KRAKEN-TRANSLATE ON OUTPUT            ########################"
echo "########################################################################################################"
echo " "

for file in ${out}/*.krak; do
	${KRAKEN}/kraken-translate --db ${KRAKEN_DB} $file > ${file}en
	if [[ ! -s ${file}en ]] ; then echo "Something went wrong with running kraken-translate... Exiting."; exit 1; fi
	#rm $file
done




echo " "
echo "########################################################################################################"
echo "########################            MAKING KRONAGRAM OF ALL FILES               ########################"
echo "########################################################################################################"
echo " "

for file in ${out}/*.kraken; do
        ${SOFT}/scripts/kraken_to_krona.py $file > ${file%.*}.krona
	if [[ ! -s ${file%.*}.krona ]]; then echo "Something went wrong with making krona file from kraken file. Exiting..."; exit 1; fi
done

#make joined krona if there are multiple input samples
ct=$(ls -1 ${out}/| grep krona | wc -l)
if [ $ct -gt 1 ]; then
	cat ${out}/*.krona >> ${out}/ALL.krona
fi



ktImportText -o ${out}/kronagram.html ${out}/*krona





echo " "
echo "########################################################################################################"
echo "########################         FINISHED RUNNING KRAKEN PIPELINE!!!            ########################"
echo "########################################################################################################"
echo " "














