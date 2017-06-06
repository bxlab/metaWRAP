#!/bin/bash -l
#SBATCH
#SBATCH --job-name=QC_pipeline
#SBATCH --partition=shared
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G

echo " "
echo "##############################################################################################################################################################"
echo "#"
echo "# INPUT: sbatch QC_PIPELINE.sh reads_1.fastq reads_2.fastq"
echo "# This script is meant to be a comprehensive solution to QC new HiSeq reads in preparation for assembly, and other operations."
echo "# The main things this pipeline accomplishes are read trimming based on quality scores, and removal of human sequences."
echo "# The script also produces a FASTQC report before and after the procedures."
echo "# Requires bmtagger to be installed in the PATH! (by conda for example)"
echo "#"
echo "# Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses." 
echo "# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu."
echo "# "
echo "##############################################################################################################################################################"
echo " "



# Don't forget to se tthe correct path to the SOFTWARE folder, containg all the needed packages:
SOFT="/home-2/guritsk1@jhu.edu/scratch/METAGENOME_PIPELINE/SOFTWARE"
# this path needs to have an indexed human genome (see bmtagger guide on their website). This insludes files hg38.bitmask and hg38.srprism.*
BMTAGGER_DB="/home-2/guritsk1@jhu.edu/scratch/Programs/bmtagger_db"
threads=24




# Checks for correct script input
if [ "$#" -ne 2 ]; then
echo "************************************************************************************************************************"
echo "*****************************************    ERROR: INVALID INPUT     **************************************************"
echo " "
echo "Please run this script like this: ./QC_PIPELINE sample_A.fastq sampleA_2.fastq"
echo "The number of files must be two, and they must be paired-end."
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





mkdir READS_QC
sample=${1%_*}
mkdir READS_QC/$sample

echo " "
echo "########################################################################################################"
echo "########################                 MAKING PRE-QC REPORT                   ########################"
echo "########################################################################################################"
echo " "

mkdir READS_QC/${sample}/pre-QC_report
${SOFT}/fastqc/fastqc -t $threads -o READS_QC/${sample}/pre-QC_report -f fastq $1 $2
rm READS_QC/${sample}/pre-QC_report/*zip
echo "####  pre-qc report saved to: READS_QC/${sample}/pre-QC_report"
echo " "
echo " "




echo " "
echo "########################################################################################################"
echo "########################                 RUNNING TRIM-GALORE                    ########################"
echo "########################################################################################################"
echo " "

${SOFT}/trim_galore --no_report_file --paired -o READS_QC/$sample $1 $2
# Fix the naming of the trimmed reads files:
for i in READS_QC/$sample/*fq; do 
	tmp1=${i%_*}; 
	tmp2=${tmp1%_*}; 
	mv $i ${tmp2}.trimmed.fastq;
	if [[ ! -s ${tmp2}.trimmed.fastq ]]; then echo "Something went wrong with trimming the reads. Exiting."; exit 1; fi
done
echo "####  Trimmed reads saved to: READS_QC/${sample}/*.trimmed.fastq"
echo " "
echo " "





echo " "
echo "########################################################################################################"
echo "########################               REMOVING HUMAN SEQUENCES                 ########################"
echo "########################################################################################################"
echo " "

mkdir READS_QC/${sample}/bmtagger_tmp
bmtagger.sh -b ${BMTAGGER_DB}/hg38.bitmask -x ${BMTAGGER_DB}/hg38.srprism -T READS_QC/${sample}/bmtagger_tmp -q1\
 -1 READS_QC/${sample}/${sample}_1.trimmed.fastq\
 -2 READS_QC/${sample}/${sample}_2.trimmed.fastq\
 -o READS_QC/${sample}/${sample}.bmtagger.list
if [[ ! -s READS_QC/${sample}/${sample}.bmtagger.list ]]; then echo "WARNING: Something went wrong with finding human reads (unlikely that there are none...). Just sayin'"; fi

rm READS_QC/${sample}/bmtagger_tmp

echo "Now sorting out found human reads from the main fastq files..."
${SOFT}/scripts/skip_human_reads.py READS_QC/${sample}/${sample}.bmtagger.list READS_QC/${sample}/${sample}_1.trimmed.fastq > READS_QC/${sample}/${sample}_1.clean.fastq
${SOFT}/scripts/skip_human_reads.py READS_QC/${sample}/${sample}.bmtagger.list READS_QC/${sample}/${sample}_2.trimmed.fastq > READS_QC/${sample}/${sample}_2.clean.fastq

echo "Now sorting out found human reads and putting them into a new file... for science..."
${SOFT}/scripts/select_human_reads.py READS_QC/${sample}/${sample}.bmtagger.list READS_QC/${sample}/${sample}_1.trimmed.fastq > READS_QC/${sample}/human_reads_1.fastq
${SOFT}/scripts/select_human_reads.py READS_QC/${sample}/${sample}.bmtagger.list READS_QC/${sample}/${sample}_2.trimmed.fastq > READS_QC/${sample}/human_reads_2.fastq
if [[ ! -s READS_QC/${sample}/${sample}_1.clean.fastq ]]; then echo "Something went wrong with removing human reads with bmtagger. Exiting."; exit 1; fi

echo "####  Contamination-free and trimmed reads are stored in: READS_QC/${sample}/${sample}_1-2.clean.fastq"
echo " "
echo " "





echo " "
echo "########################################################################################################"
echo "########################                 MAKING POST-QC REPORT                  ########################"
echo "########################################################################################################"
echo " "

mkdir READS_QC/${sample}/post-QC_report
${SOFT}/fastqc/fastqc -t $threads -o READS_QC/${sample}/post-QC_report -f fastq READS_QC/${sample}/${sample}_1.clean.fastq READS_QC/${sample}/${sample}_2.clean.fastq
rm READS_QC/${sample}/post-QC_report/*zip
echo "####  pre-qc report saved to: READS_QC/${sample}/post-QC_report"
echo " "
echo " "


# Final copy of the good reads
mv READS_QC/${sample}/${sample}_1.clean.fastq ${sample}_1.clean.fastq
mv READS_QC/${sample}/${sample}_2.clean.fastq ${sample}_2.clean.fastq

# Remove intermediate files
#rm READS_QC/${sample}/*fastq


echo " "
echo "########################################################################################################"
echo "########################              READ QC PIPELINE COMPLETE!!!              ########################"
echo "########################################################################################################"
echo " "








