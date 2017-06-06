#!/bin/bash -l
#SBATCH
#SBATCH --job-name=ASSEMBLY_pipeline
#SBATCH --partition=lrgmem
#SBATCH --time=160:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=900G

echo " "
echo "##############################################################################################################################################################"
echo "#"
echo "# INPUT: sbatch ASSEMBLY_PIPELINE.sh assembly_output sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq ..."
echo "# This script is meant to be a comprehensive solution for producing the best metagenomic assembly given paired end reads from one or more samples."
echo "# Ideally it should take in fully QC'd reads. First, the reads are assembled with metaSPAdes3.10, then all the reads that did not map back to the"
echo "# contigs are re-assembled with MEGAHIT (which works better on lower coverage contigs. The resulting assemblies are combined, sorted, and short "
echo "# contigs are removed. The finall assembly is then QCed by QUAST. (Please indicate path to quast at the sart of script)"
echo "#"
echo "# Author of pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses."
echo "# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu."
echo "# "
echo "##############################################################################################################################################################"
echo " "




# Don't forget to set the correct path to the SOFTWARE folder, containg all the needed packages. This must be an absolute path:"
SOFT="/home-2/guritsk1@jhu.edu/scratch/METAGENOME_PIPELINE/SOFTWARE"
QUAST="/home-2/guritsk1@jhu.edu/scratch/Programs/quast-4.0"
threads=48
mem=850




# Checks for correct script input
if [ "$#" -lt 3 ]; then
echo "************************************************************************************************************************"
echo "*****************************************    ERROR: INVALID INPUT     **************************************************"
echo " "
echo "Please run this script like this: ./ASSEMBLY_PIPELINE output_folder sample_A.fastq sampleA_2.fastq [sampleX_1 sampleX_2]"
echo "The number of files must be odd, since the reads must be paired and labeled *_1.fastq or *_2.fastq."
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








echo " "
echo "########################################################################################################"
echo "########################            COMBINING SAMPLES (IF MULTIPLE)             ########################"
echo "########################################################################################################"
echo " "

out=$1
mkdir $out

# If there are several pairs of reads passed, they are combined into two files of F and R reads
if [ "$#" -ne 3 ]; then
for num in "$@"; do 
if [[ $num == *"_1."*"fastq"* ]]; then cat $num >> ${out}/joined_reads_1.fastq; fi
if [[ $num == *"_2."*"fastq"* ]]; then cat $num >> ${out}/joined_reads_2.fastq; fi
done
reads_1=${out}/joined_reads_1.fastq
reads_2=${out}/joined_reads_2.fastq
fi

# If there are exactly two files passed, those will be the paired reads
if [ "$#" -eq 3 ]; then
reads_1=$2
reads_2=$3
fi

echo " "
echo "### Using reads $reads_1 and $reads_2 for assembly. ####"
echo " "

if [ ! -f "$reads_1" ] || [ ! -f "$reads_2" ]; then echo "Something went wrong with serring up the reads for assembly. Exiting."; exit 1; fi



echo " "
echo "########################################################################################################"
echo "########################               ASSEMBLING WITH METASPADES               ########################"
echo "########################################################################################################"
echo " "

${SOFT}/SPAdes-3.10.1-Linux/bin/metaspades.py --tmp-dir ${out}/metaspades.tmp -t $threads -m $mem -o ${out}/metaspades -1 $reads_1 -2 $reads_2 
rm ${out}/metaspades.tmp
if [ ! -f "${out}/metaspades/scaffolds.fasta" ]; then echo "Something went wrong with metaSPAdes assembly. Exiting."; exit 1; fi



echo " "
echo "########################################################################################################"
echo "########################              SORTING OUT UNASSEMBLED READS             ########################"
echo "########################################################################################################"
echo " "

#take only the scaffolds over 1.5kb (the minimum needed for binning with metaBAT):
${SOFT}/scripts/rm_short_contigs.py 1500 ${out}/metaspades/scaffolds.fasta > ${out}/metaspades/long_scaffolds.fasta
${SOFT}/bwa index ${out}/metaspades/long_scaffolds.fasta

#sort out and store reads that dont map back to the assembly:
${SOFT}/bwa mem -t $threads ${out}/metaspades/long_scaffolds.fasta $reads_1 $reads_2 | grep -v NM:i: > ${out}/unused_by_metaspades.sam
${SOFT}/scripts/sam_to_fastq.py ${out}/unused_by_metaspades.sam > ${out}/unused_by_metaspades.fastq
rm ${out}/unused_by_metaspades.sam

if [[ -s ${out}/unused_by_metaspades.sam ]]; then echo "Something went wrong with pulling out unassembled reads. Exiting."; exit 1; fi



echo " "
echo "########################################################################################################"
echo "########################        ASSEMBLING REMAINDER READS WITH MEGAHIT         ########################"
echo "########################################################################################################"
echo " "

${SOFT}/megahit/megahit -r ${out}/unused_by_metaspades.fastq -o ${out}/megahit -t $threads -m ${mem}000000000
if [ ! -f "${out}/megahit/final.contigs.fa" ]; then echo "Something went wrong with reassembling with Megahit. Exiting."; exit 1; fi



echo " "
echo "########################################################################################################"
echo "########################         COMBINE AND FORMAT THE TWO ASSEMBLIES          ########################"
echo "########################################################################################################"
echo " "

${SOFT}/scripts/fix_megahit_contig_naming.py ${out}/megahit/final.contigs.fa > ${out}/megahit/fixed.contigs.fa
${SOFT}/scripts/rm_short_contigs.py 500 ${out}/megahit/fixed.contigs.fa > ${out}/megahit/long.contigs.fa
cp ${out}/metaspades/long_scaffolds.fasta ${out}/combined_assembly.fasta
cat ${out}/megahit/long.contigs.fa >> ${out}/combined_assembly.fasta
${SOFT}/scripts/sort_contigs.py ${out}/combined_assembly.fasta > ${out}/final_assembly.fasta
rm ${out}/combined_assembly.fasta
if [[ ! -s ${out}/${out}/combined_assembly.fasta ]]; then echo "Something went wrong with joining the assemblies. Exiting."; exit 1; fi



echo " "
echo "########################################################################################################"
echo "########################             RUNNING ASSEMBLY QC WITH QUAST             ########################"
echo "########################################################################################################"
echo " "

cp ${out}/metaspades/scaffolds.fasta ${out}/metaspades/metaspades_assembly.fasta
${QUAST}/quast.py -t $threads -o ${out}/QUAST_out -m 500 ${out}/final_assembly.fasta ${out}/metaspades/metaspades_assembly.fasta




echo " "
echo "########################################################################################################"
echo "########################             ASSEMBLY PIPELINE COMPLETED!!!             ########################"
echo "########################################################################################################"
echo " "





























