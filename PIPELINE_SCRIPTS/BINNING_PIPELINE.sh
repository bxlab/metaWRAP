#!/bin/bash -l
#SBATCH
#SBATCH --job-name=BINNING_pipeline
#SBATCH --partition=lrgmem
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=900G

echo " "
echo "##############################################################################################################################################################"
echo "#"
echo "# INPUT: sbatch BINNING_PIPELINE.sh assembly.fa sampleA_1.fastq sampleA_2.fastq [ sampleX_1.fastq sampleX_2.fastq ]"
echo "# This script is meant to be run on the outputs of the QC_PIPELINE and ASSSEMBLY_PIPELINE to split the assembly contigs into metagenomic bins."
echo "# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly."
echo "# The more samples, the better the binning. "
echo "# "
echo "# The script uses metaBAT to bin the contigs, then uses bwa to recruit reads back to the bins and ressembles them with SPAdes. It then uses KRAKEN "
echo "# to assign taxonomy to each bin and producing a joint kronagram of all the bins. Finally, it uses CheckM to test the contamination and completion"
echo "# of the bins, and sepperates bins with \>20% completion and \<10% contamintation. "
echo " "
echo "# NOTE: Samtools, CheckM, and KRAKEN require instalation, and be sure the configure the right path to the KRAKEN folder and the KRAKEN database."
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
CHECKM_PATH="/home-2/guritsk1@jhu.edu/bin"
threads=48
mem=850



# Checks for correct script input
if [ "$#" -lt 3 ]; then
echo "************************************************************************************************************************"
echo "*****************************************    ERROR: INVALID INPUT     **************************************************"
echo " "
echo "Please run this script like this: ./BINNING_PIPELINE assembly.fa sample_A.fastq sampleA_2.fastq [sampleX_1 sampleX_2]"
echo "The number of read files must be even, since the reads must be paired"
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


#Checks for present CheckM installation
if [ ! -f "${CHECKM_PATH}/checkm" ]; then
echo "************************************************************************************************************************"
echo "*******************************           ERROR: CHECKM NOT FOUND!           *******************************************"
echo " "
echo "Depending on your system/cluster, CheckM may be a pain to install... Good luck! Its worth it. "
echo " "
echo "************************************************************************************************************************"
exit 1
fi







echo " "
echo "########################################################################################################"
echo "########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################"
echo "########################################################################################################"
echo " "

# setting up the output folder
sample=${1##*/}
binning=${sample%.*}.binning
mkdir $binning
mkdir ${binning}/alignments
cp $1 ${binning}/alignments/assembly.fa





# Index the assembly
${SOFT}/bwa index ${binning}/alignments/assembly.fa
if [[ ! -s ${binning}/alignments/assembly.fa.amb ]] ; then echo "Something went wrong with indexing the assembly.\
 ${binning}/alignments/assembly.fa.amb does not exits! Exiting."; exit 1; fi

# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	if [[ $num == *"_1."*"fastq"* ]]; then 
		reads_1=$num
		reads_2="${reads_1/_1/_2}"
		tmp=${reads_1##*/}
		sample=${tmp%_*}

		echo "####       Processing $reads_1 and $reads_2         ####"
		${SOFT}/bwa mem -t $threads ${binning}/alignments/assembly.fa $reads_1 $reads_2 | samtools view -bS\
		 | samtools sort -O BAM --threads $threads -o ${binning}/alignments/${sample}.bam -
		if [[ ! -s ${binning}/alignments/${sample}.bam ]]; then echo "Something went wrong with aligning reads to contigs. Exiting."; exit 1; fi
	fi
done


echo " "
echo "########################################################################################################"
echo "########################                   RUNNING METABAT2                     ########################"
echo "########################################################################################################"
echo " "

${SOFT}/metabat2/jgi_summarize_bam_contig_depths --outputDepth ${binning}/depth.txt ${binning}/alignments/*.bam
${SOFT}/metabat2/metabat2 -i ${binning}/alignments/assembly.fa -a ${binning}/depth.txt -o ${binning}/metabat2_bins/bin -m 1500 --unbinned
if [[ ! -s ${binning}/metabat2_bins/bin.unbinned.fa ]]; then echo "Something went wrong with running MetaBAT. Exiting"; exit 1; fi
mv ${binning}/depth.txt ${binning}/alignments/



echo " "
echo "########################################################################################################"
echo "########################        RECRUITING READS TO BINS FOR REASSEMBLY         ########################"
echo "########################################################################################################"
echo " "

mkdir ${binning}/reassembly

echo "Combining reads..."
for num in "$@"; do
	if [[ $num == *"_1."*"fastq"* ]]; then
		reads_1=$num
		reads_2="${reads_1/_1/_2}"
		cat $reads_1 >> ${binning}/reassembly/joined_reads_1.fastq
		cat $reads_2 >> ${binning}/reassembly/joined_reads_2.fastq
	fi
done
if [[ ! -s ${binning}/reassembly/joined_reads_2.fastq ]] ; then echo "Something went wrong with combining reads for reassembly. Exiting."; exit 1; fi


echo "Move over bin files, fix their naming, and index them in preperation for alignment"
if [ -d ${binning}/reassembly/bins ]; then rm -r ${binning}/reassembly/bins; fi
mkdir ${binning}/reassembly/bins
cp ${binning}/metabat2_bins/* ${binning}/reassembly/bins/
rm ${binning}/reassembly/bins/*unbinned*
echo "Making sure that the extension is .fa"
for i in ${binning}/reassembly/bins/*; do mv $i ${i%.*}.fa; done
for i in ${binning}/reassembly/bins/*.fa; do ${SOFT}/bwa index $i; done



echo "Align reads back to bins, and save aligned reads"
mkdir ${binning}/reassembly/reads
for i in ${binning}/reassembly/bins/*.fa; do 
	base=${i##*/};
	echo "#### ALIGNING READS ONTO ${i##*/} FILE #####"
	${SOFT}/bwa mem -t $threads $i ${binning}/reassembly/joined_reads_1.fastq ${binning}/reassembly/joined_reads_2.fastq\
	| grep NM:i: | cut -f1 | uniq > ${binning}/reassembly/reads/${base%.*}.list;

	echo "#### EXTRACTING PAIRED-END READS INTO NEW FILE TO REASSEMBLE ${i##*/} #####"
	# save the aligned sequences in new files
	${SOFT}/scripts/filter_out_fastq_reads.py ${binning}/reassembly/reads/${base%.*}.list ${binning}/reassembly/joined_reads_1.fastq\
	 > ${binning}/reassembly/reads/${base%.*}_1.fastq
	${SOFT}/scripts/filter_out_fastq_reads.py ${binning}/reassembly/reads/${base%.*}.list ${binning}/reassembly/joined_reads_2.fastq\
	 > ${binning}/reassembly/reads/${base%.*}_2.fastq

	if [[ ! -s ${binning}/reassembly/reads/${base%.*}_2.fastq ]]; then echo "Something went wrong with extracting reads for reassembling the bins. Exiting."; exit 1; fi
done






echo " "
echo "########################################################################################################"
echo "########################             REASSEMBLING BINS WITH SPADES              ########################"
echo "########################################################################################################"
echo " "

mkdir ${binning}/reassembly/reassemblies
for bin in ${binning}/reassembly/reads/*.list; do
	bin_name=${bin##*/}
	echo "#### NOW PROCESSING ${bin_name%.*} ####"
	${SOFT}/SPAdes-3.10.1-Linux/bin/spades.py -t $threads -m $mem --tmp /tmp --careful \
	--trusted-contigs ${binning}/reassembly/bins/${bin_name%.*}.fa \
	-1 ${binning}/reassembly/reads/${bin_name%.*}_1.fastq \
	-2 ${binning}/reassembly/reads/${bin_name%.*}_2.fastq \
	-o ${binning}/reassembly/reassemblies/${bin_name%.*}
	if [[ ! -s ${binning}/reassembly/reassemblies/${bin_name%.*}/scaffolds.fasta ]]; then echo "Something went wrong with reassembling ${bin_name%.*}. Exiting."; exit 1; fi
done













echo " "
echo "########################################################################################################"
echo "########################    ALIGNING READS TO PREPARE TO RE-BIN ASSEMBLIES      ########################"
echo "########################################################################################################"
echo " "

mkdir ${binning}/rebinning
for bin in ${binning}/reassembly/reassemblies/*; do
	bin_name=${bin##*/}
	echo "#### NOW PROCESSING $bin_name ####"

	mkdir ${binning}/rebinning/${bin_name}
	cp ${bin}/scaffolds.fasta ${binning}/rebinning/${bin_name}/reassembled_bin.fa
	# Index the assembly
	${SOFT}/bwa index ${binning}/rebinning/${bin_name}/reassembled_bin.fa

	# If there are several pairs of reads passed, they are processed sepperately
	for num in "$@"; do
	        if [[ $num == *"_1."*"fastq"* ]]; then
	                reads_1=$num
	                reads_2="${reads_1/_1/_2}"
	                tmp=${reads_1##*/}
	                sample=${tmp%_*}
	
	                echo "####       Processing $reads_1 and $reads_2         ####"
	                ${SOFT}/bwa mem -t $threads ${binning}/rebinning/${bin_name}/reassembled_bin.fa $reads_1 $reads_2 |\
 			samtools view -bS - | samtools sort -O BAM --threads $threads -o ${binning}/rebinning/${bin_name}/${sample}.bam -
			if [[ ! -s ${binning}/rebinning/${bin_name}/${sample}.bam ]]; then echo "Something went wrong with aligning\
			 reads to reassembled bins to prepare for re-binning. Exiting."; exit 1; fi
	        fi
	done
done






echo " "
echo "########################################################################################################"
echo "########################       REMOVING ERRONEOUS CONTIGS WITH METABAT          ########################"
echo "########################################################################################################"
echo " "

for bin in ${binning}/rebinning/*; do
	bin_name=${bin##*/}
	echo "#### NOW PROCESSING $bin_name ####"
	${SOFT}/metabat/jgi_summarize_bam_contig_depths --outputDepth ${binning}/rebinning/${bin_name}/depth.txt ${binning}/rebinning/${bin_name}/*.bam
	echo "#### NOW BINNING $bin_name ####"
	${SOFT}/metabat2/metabat2 -m 1500 --maxP 99 --minS 40 --maxEdges 300 \
	 -i ${binning}/rebinning/${bin_name}/reassembled_bin.fa\
	 -a ${binning}/rebinning/${bin_name}/depth.txt\
	 -o ${binning}/rebinning/${bin_name}/metabat/bin
	if [[ ! -s ${binning}/rebinning/${bin_name}/metabat/bin.1.fa ]]; then echo "WARNING: Something went wrong with re-binning ${bin_name} \
	 Either the pipeline is failing or the bin was lost in the reassembly process"; fi
done

# copy over the best bin from each re-binning
mkdir ${binning}/reassembled_bins
for bin in ${binning}/rebinning/*; do
	best_bin=$(eval du -a ${bin}/metabat | sort -n -r | head -2 | tail -1 | cut -f2)
	#best_bin=$(eval du -a ${bin}/metabat | grep -v unbinned | sort -n -r | head -2 | tail -1 | cut -f2)
	bin_base=${bin##*/}
	echo "cp $best_bin ${binning}/reassembled_bins/${bin_base}.fa"
	cp $best_bin ${binning}/reassembled_bins/${bin_base}.fa


done







echo " "
echo "########################################################################################################"
echo "########################              RUNNING KRAKEN ON FINAL BINS              ########################"
echo "########################################################################################################"
echo " "

mkdir ${binning}/KRAKEN
for bin in ${binning}/reassembled_bins/*.fa; do
        base_name=${bin##*/}
        ${KRAKEN}/kraken --db ${KRAKEN_DB} --fasta-input --threads $threads\
         --output ${binning}/KRAKEN/${base_name%.*}.krak $bin
        if [[ ! -s ${binning}/KRAKEN/${base_name%.*}.krak ]] ; then echo "WARNING: Something went wrong with running kraken on $bin ..."; fi

        ${KRAKEN}/kraken-translate --db ${KRAKEN_DB} ${binning}/KRAKEN/${base_name%.*}.krak > ${binning}/KRAKEN/${base_name%.*}.kraken
        if [[ ! -s ${binning}/KRAKEN/${base_name%.*}.kraken ]] ; then echo "WARNING: Something went wrong with running kraken-translate on $bin"; fi
        rm ${binning}/KRAKEN/${base_name%.*}.krak
done


echo " "
echo "########################################################################################################"
echo "########################            MAKING KRONAGRAM OF ALL FILES               ########################"
echo "########################################################################################################"
echo " "

for file in ${binning}/KRAKEN/*.kraken; do
        ${SOFT}/scripts/kraken_to_krona.py $file > ${file%.*}.krona
	if [[ ! -s ${file%.*}.krona ]]; then echo "WARNING: Something went wrong with making krona file from kraken file $file ..."; fi
done


ktImportText -o ${binning}/KRAKEN/kronagram.html ${binning}/KRAKEN/*krona







echo " "
echo "########################################################################################################"
echo "########################             RUN CHECKM ON RESULTING BINS               ########################"
echo "########################################################################################################"
echo " "

#NOTE: The following code is very system dependant. I only have to do this becasue I am running this on MARCC:
module load python/2.7.10
#this just loaded the correct python version to run my compilation of checkm...

rm -r ${binning}/metabat2_bins.checkm ${binning}/reassembled_bins.checkm

echo "Running CheckM on original metaBAT2 bins"
python ${CHECKM_PATH}/checkm lineage_wf -x fa ${binning}/metabat2_bins ${binning}/metabat2_bins.checkm -t $threads
if [[ ! -s ${binning}/metabat2_bins.checkm/storage/bin_stats_ext.tsv ]]; then echo "Something went wrong with running CheckM. Exiting..."; exit 1; fi
echo "Making CheckM plot of original metaBAT bins"
python ${CHECKM_PATH}/checkm bin_qa_plot -x fa ${binning}/metabat2_bins.checkm ${binning}/metabat2_bins ${binning}/metabat2_bins.plot
if [[ ! -s ${binning}/metabat2_bins.plot/bin_qa_plot.png ]]; then echo "Something went wrong with making the CheckM plot. Exiting."; exit 1; fi



echo "Running CheckM on reassembled metaBAT bins"
python ${CHECKM_PATH}/checkm lineage_wf -x fa ${binning}/reassembled_bins ${binning}/reassembled_bins.checkm -t $threads
if [[ ! -s ${binning}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then echo "Something went wrong with running CheckM. Exiting..."; exit 1; fi
echo "Making CheckM plot of reassembled metaBAT bins"
python ${CHECKM_PATH}/checkm bin_qa_plot -x fa ${binning}/reassembled_bins.checkm ${binning}/reassembled_bins ${binning}/reassembled_bins.plot
if [[ ! -s ${binning}/reassembled_bins.plot/bin_qa_plot.png ]]; then echo "Something went wrong with making the CheckM plot. Exiting."; exit 1; fi



echo "Finalizing..."
${SOFT}/scripts/summarize_checkm.py ${binning}/metabat2_bins.checkm/storage/bin_stats_ext.tsv > ${binning}/metabat2_bins.stats
${SOFT}/scripts/summarize_checkm.py ${binning}/reassembled_bins.checkm/storage/bin_stats_ext.tsv > ${binning}/reassembled_bins.stats
mv ${binning}/metabat2_bins.plot/bin_qa_plot.png ${binning}/metabat2_bins.png
mv ${binning}/reassembled_bins.plot/bin_qa_plot.png ${binning}/reassembled_bins.png
rm -r ${binning}/reassembled_bins.plot ${binning}/metabat2_bins.plot


echo " "
echo "########################################################################################################"
echo "########################            BINNING PIPELINE FINISHED!!!                ########################"
echo "########################################################################################################"
echo " "


























