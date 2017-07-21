#!/bin/bash -l

##############################################################################################################################################################
#
# This script is meant to be run on the outputs of assembly.sh pipeline to split the assembly contigs into metagenomic bins.
# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.
# The more samples, the better the binning. 
#
# The script uses metaBAT to bin the contigs, then uses bwa to recruit reads back to the bins and ressembles them with SPAdes. It then uses KRAKEN 
# to assign taxonomy to each bin and producing a joint kronagram of all the bins. Finally, it uses CheckM to test the contamination and completion
# of the bins, and sepperates bins with \>20% completion and \<10% contamintation. 
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################


help_message () {
	echo "Usage: ./binning.sh [options] -a assembly.fa -o output_dir readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]"
	echo "Options:"
	echo ""
	echo "	-a STR          metagenomic assembly"
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT          memory in GB (default=4)"
	echo "";}
# function to print out error messages
error () {
	echo ""; echo "************************************************************************************************************************"
	echo "*****************************************            ERROR!           **************************************************"
	echo $1
	echo "************************************************************************************************************************"; echo ""; exit 1; }
# function to print out warning messages
warning () {
	echo ""; echo  "************************************************************************************************************************"
	echo "*****************************************           WARNING!          **************************************************"
	echo $1
	echo "************************************************************************************************************************"; echo ""; }
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

threads=1; mem=4; out="false"; ASSEMBLY="false"; 
# Load in options
while getopts ht:m:o:a: option; do
	case "${option}" in
		h) help_message; exit 1;;
		t) threads=${OPTARG};;
		m) mem=${OPTARG};;
		o) out=${OPTARG};;
		a) ASSEMBLY=$OPTARG;;
	esac
done




########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$#" -lt 2 ] || [ "$ASSEMBLY" = "false" ] ; then 
	help_message; exit 1
fi

# check for at least one pair of read fastq files:
F="no"; R="no"
for num in "$@"; do
	if [[ $num == *"_1.fastq" ]]; then F="yes"; fi
	if [[ $num == *"_2.fastq" ]]; then R="yes"; fi
done
if [ $F = "no" ] || [ $R = "no" ]; then
	comm "Unable to find proper fastq read pair in the format *_1.fastq and *_2.fastq"
	help_message; exit 1
fi

#determine number of fastq read files provided:
num_of_F_read_files=$(for I in "$@"; do echo $I | grep _1.fastq; done | wc -l)
num_of_R_read_files=$(for I in "$@"; do echo $I | grep _2.fastq; done | wc -l)

comm "$num_of_F_read_files forward and $num_of_R_read_files reverse read files detected"
if [ ! $num_of_F_read_files == $num_of_R_read_files ]; then error "Number of F and R reads must be the same!"; fi


# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi

# Checks for KRAKEN database 
if [ ! -d "$KRAKEN_DB" ]; then
	error "The folder $KRAKEN_DB doesnt exist. Please look and the the script and manually set the correct path to the KRAKEN standard database. \
	If you do not have it yet, you can download it like this: ${KRAKEN}/kraken-build --standard --threads $threads --db ${KRAKEN}/kraken/KRAKEN_DB"
fi

########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################




echo " "
echo "########################################################################################################"
echo "########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################"
echo "########################################################################################################"
echo " "

# setting up the output folder
mkdir $out ${out}/work_files ${out}/work_files/alignments
comm "making copy of assembly file $ASSEMBLY"
cp $ASSEMBLY ${out}/work_files/alignments/assembly.fa
tmp=${ASSEMBLY##*/}
sample=${tmp%.*}



# Index the assembly
comm "Indexing assembly file"
bwa index ${out}/work_files/alignments/assembly.fa
if [[ ! -s ${out}/work_files/alignments/assembly.fa.amb ]] ; then error "Something went wrong with indexing the assembly.\
 ${out}/alignments/assembly.fa.amb does not exits! Exiting."; fi


# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	if [[ $num == *"_1.fastq"* ]]; then 
		reads_1=$num
		reads_2=${num%_*}_2.fastq
		if [ "$reads_2" = "false" ]; then error "$reads_2 does not exist. Exiting..."; fi

		tmp=${reads_1##*/}
		sample=${tmp%_*}

		comm "Aligning $reads_1 and $reads_2 back to assembly"
		bwa mem -t $threads ${out}/work_files/alignments/assembly.fa $reads_1 $reads_2 | samtools view -bS\
		 | samtools sort -O BAM --threads $threads -o ${out}/work_files/alignments/${sample}.bam -
		if [[ ! -s ${out}/work_files/alignments/${sample}.bam ]]; then error "Something went wrong with aligning reads to contigs. Exiting."; fi
	fi
done


echo " "
echo "########################################################################################################"
echo "########################                   RUNNING METABAT2                     ########################"
echo "########################################################################################################"
echo " "

jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/depth.txt ${out}/work_files/alignments/*.bam
metabat2 -i ${out}/work_files/alignments/assembly.fa -a ${out}/work_files/depth.txt\
 -o ${out}/metabat2_bins/bin -m 1500 --unbinned

if [[ ! -s ${out}/metabat2_bins/bin.unbinned.fa ]]; then error "Something went wrong with running MetaBAT. Exiting"; fi



echo " "
echo "########################################################################################################"
echo "########################             RUNNING KRAKEN ON ORIGINAL BINS            ########################"
echo "########################################################################################################"
echo " "

if False; then

mkdir ${out}/kraken
for bin in ${out}/metabat2_bins/*.fa; do
        base_name=${bin##*/}
	comm "Running kraken on $base_name"
        kraken --db ${KRAKEN_DB} --fasta-input --threads $threads\
         --output ${out}/kraken/${base_name%.*}.krak $bin
        if [[ ! -s ${out}/kraken/${base_name%.*}.krak ]] ; then error "Something went wrong with running kraken on $bin ..."; fi

        kraken-translate --db ${KRAKEN_DB} ${out}/kraken/${base_name%.*}.krak > ${out}/kraken/${base_name%.*}.kraken
        if [[ ! -s ${out}/kraken/${base_name%.*}.kraken ]] ; then error "Something went wrong with running kraken-translate on $bin"; fi
        rm ${out}/kraken/${base_name%.*}.krak
done






echo " "
echo "########################################################################################################"
echo "########################            MAKING KRONAGRAM OF ALL FILES               ########################"
echo "########################################################################################################"
echo " "

for file in ${out}/kraken/*.kraken; do
        ${SOFT}/scripts/kraken_to_krona.py $file > ${file%.*}.krona
	if [[ ! -s ${file%.*}.krona ]]; then error "Something went wrong with making krona file from kraken file $file ..."; fi
done

ktImportText -o ${out}/kraken/${out%.*}_kraken.html ${out}/kraken/*krona


fi


if false; then

echo " "
echo "########################################################################################################"
echo "########################             RUN CHECKM ON RESULTING BINS               ########################"
echo "########################################################################################################"
echo " "

if [[ -d ${out}/metabat2_bins.checkm ]]; then rm -r ${out}/metabat2_bins.checkm; fi

comm "Running CheckM on original metaBAT2 bins"
checkm lineage_wf -x fa ${out}/metabat2_bins ${out}/metabat2_bins.checkm -t $threads
if [[ ! -s ${out}/metabat2_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
comm "Making CheckM plot of original metaBAT bins"
checkm bin_qa_plot -x fa ${out}/metabat2_bins.checkm ${out}/metabat2_bins ${out}/metabat2_bins.plot
if [[ ! -s ${out}/metabat2_bins.plot/bin_qa_plot.png ]]; then error "Something went wrong with making the CheckM plot. Exiting."; fi


comm "Finalizing..."
${SOFT}/summarize_checkm.py ${out}/metabat2_bins.checkm/storage/bin_stats_ext.tsv > ${out}/metabat2_bins.stats
mv ${out}/metabat2_bins.plot/bin_qa_plot.png ${out}/metabat2_bins.png
rm -r ${out}/metabat2_bins.plot
mv ${out}/metabat2_bins.checkm ${out}/work_files


fi


echo " "
echo "########################################################################################################"
echo "########################        RECRUITING READS TO BINS FOR REASSEMBLY         ########################"
echo "########################################################################################################"
echo " "

mkdir ${out}/work_files/reassembly

if [ $num_of_F_read_files gt 1 ]; then
	comm "Combining reads..."
	for num in "$@"; do
		if [[ $num == *"_1.fastq"* ]]; then
			reads_1=$num
			reads_2=${reads_1%_*}_2.fastq
			cat $reads_1 >> ${out}/work_files/reassembly/joined_reads_1.fastq
			cat $reads_2 >> ${out}/work_files/reassembly/joined_reads_2.fastq
		fi
	done
	
	F_reads=${out}/work_files/reassembly/joined_reads_1.fastq
        R_reads=${out}/work_files/reassembly/joined_reads_2.fastq
	
	if [[ ! -s $F_reads ]] ; then
		error "Something went wrong with combining reads for reassembly. Exiting." 
	fi

else
	for file in "$@"; do
                if [[ $file == *"_1.fastq" ]]; then F_reads=$file; fi
		if [[ $file == *"_2.fastq" ]]; then R_reads=$file; fi
	done
fi



comm "Moving over bin files, fixing their naming, and indexing them in preperation for alignment"
rm -r ${out}/work_files/reassembly/bins
mkdir ${out}/work_files/reassembly/bins
cp ${out}/metabat2_bins/* ${out}/work_files/reassembly/bins/
comm "Making sure that the extension is .fa"
for i in ${out}/work_files/reassembly/bins/*; do mv $i ${i%.*}.fa; done
for i in ${out}/work_files/reassembly/bins/*.fa; do bwa index $i; done


comm "Align reads back to bins, and save aligned reads"
rm -r ${out}/work_files/reassembly/reads
mkdir ${out}/work_files/reassembly/reads
for i in ${out}/work_files/reassembly/bins/*.fa; do 
	base=${i##*/};
	
	comm "ALIGNING READS ONTO ${i##*/} FILE"
	#find perfectly aligning reads
	echo "Finding perfect alignments"
	bwa mem -t $threads $i $F_reads $R_reads\
	| grep NM:i:0 | cut -f1 | uniq > ${out}/work_files/reassembly/reads/${base%.*}.strict.list;

	#find all aligning reads
	echo "Finding all alignment"
	bwa mem -t $threads $i $F_reads $R_reads\
        | grep NM:i: | cut -f1 | uniq > ${out}/work_files/reassembly/reads/${base%.*}.permissive.list;


	comm "EXTRACTING PAIRED-END READS INTO NEW FILE TO REASSEMBLE ${i##*/}"
	# save the perfectly aligned sequences in new files
	echo "Extracting perfectly aligned reads for STRICT reassembly"
	${SOFT}/filter_out_fastq_reads.py ${out}/work_files/reassembly/reads/${base%.*}.strict.list $F_reads\
	 > ${out}/work_files/reassembly/reads/${base%.*}_1.strict.fastq
	${SOFT}/filter_out_fastq_reads.py ${out}/work_files/reassembly/reads/${base%.*}.strict.list $R_reads\
	 > ${out}/work_files/reassembly/reads/${base%.*}_2.strict.fastq


	# save all aligned sequences in new files
	echo "Extracting all aligned reads for PERMISSIVE reassembly"
        ${SOFT}/filter_out_fastq_reads.py ${out}/work_files/reassembly/reads/${base%.*}.permissive.list $F_reads\
         > ${out}/work_files/reassembly/reads/${base%.*}_1.permissive.fastq
        ${SOFT}/filter_out_fastq_reads.py ${out}/work_files/reassembly/reads/${base%.*}.permissive.list $R_reads\
         > ${out}/work_files/reassembly/reads/${base%.*}_2.permissive.fastq


	if [[ ! -s ${out}/work_files/reassembly/reads/${base%.*}_2.strict.fastq ]]; then
		error "Something went wrong with extracting reads for reassembling the bins. Exiting."
	fi
done




echo " "
echo "########################################################################################################"
echo "########################             REASSEMBLING BINS WITH SPADES              ########################"
echo "########################################################################################################"
echo " "

mkdir ${out}/work_files/reassembly/reassemblies
mkdir ${out}/reassembled_bins
for bin in ${out}/work_files/reassembly/reads/*_1.strict.fastq; do
	bin_name=${bin##*/}
	comm "NOW REASSEMBLING STRICT ${bin_name%_*}"
	spades.py -t $threads -m $mem --tmp /tmp --careful \
	--untrusted-contigs ${out}/work_files/reassembly/bins/${bin_name%_*}.fa \
	-1 ${out}/work_files/reassembly/reads/${bin_name%_*}_1.strict.fastq \
	-2 ${out}/work_files/reassembly/reads/${bin_name%_*}_2.strict.fastq \
	-o ${out}/work_files/reassembly/reassemblies/${bin_name%_*}.strict
	
	if [[ ! -s ${out}/work_files/reassembly/reassemblies/${bin_name%_*}.strict/scaffolds.fasta ]]; then
                warning "Something went wrong with reassembling ${bin_name%_*}.strict"
	else 
		comm "${bin_name%_*}.strict was reassembled successfully!"
	fi


	comm "NOW REASSEMBLING PERMISSIVE ${bin_name%_*}"
        spades.py -t $threads -m $mem --tmp /tmp --careful \
        --untrusted-contigs ${out}/work_files/reassembly/bins/${bin_name%_*}.fa \
        -1 ${out}/work_files/reassembly/reads/${bin_name%_*}_1.permissive.fastq \
        -2 ${out}/work_files/reassembly/reads/${bin_name%_*}_2.permissive.fastq \
        -o ${out}/work_files/reassembly/reassemblies/${bin_name%_*}.permissive
	
	
	if [[ ! -s ${out}/work_files/reassembly/reassemblies/${bin_name%_*}.permissive/scaffolds.fasta ]]; then 
		warning "Something went wrong with reassembling ${bin_name%_*}.permissive"
	else
		comm "${bin_name%_*}.permissive was reassembled successfully!"
	fi
	
done


# removing short contigs and placing reassemblies in the final folder
comm "Finalizing reassemblies"
for spades_folder in ${out}/work_files/reassembly/reassemblies/*; do
	bin_name=${spades_folder##*/}
	
	#remove shortest contigs (probably artifacts...)
	${SOFT}/rm_short_contigs.py 500\
	 ${out}/work_files/reassembly/reassemblies/${bin_name}/scaffolds.fasta\
	 > ${out}/work_files/reassembly/reassemblies/${bin_name}/long_scaffolds.fasta
	
	cp ${out}/work_files/reassembly/reassemblies/${bin_name}/long_scaffolds.fasta\
	 ${out}/reassembled_bins/${bin_name}.fa
done

if [[ ! -s ${out}/reassembled_bins/bin.unbinned.permissive.fa ]]; then
	error "Something went wrong with processing the reassembled bins and placing them into the reassembled_bins folder. Exiting." 
fi






echo " "
echo "########################################################################################################"
echo "########################             RUN CHECKM ON REASSEMBLED BINS             ########################"
echo "########################################################################################################"
echo " "

rm -r ${out}/metabat2_bins.checkm ${out}/reassembled_bins.checkm

comm "Running CheckM on reassembled metaBAT bins"
checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads
if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
comm "Making CheckM plot of reassembled metaBAT bins"
checkm bin_qa_plot -x fa ${out}/reassembled_bins.checkm ${out}/reassembled_bins ${out}/reassembled_bins.plot
if [[ ! -s ${out}/reassembled_bins.plot/bin_qa_plot.png ]]; then error "Something went wrong with making the CheckM plot. Exiting."; fi



comm "Finalizing CheckM stats and plots..."
${SOFT}/scripts/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv > ${out}/reassembled_bins.stats
mv ${out}/reassembled_bins.plot/bin_qa_plot.png ${out}/reassembled_bins.png
rm -r ${out}/reassembled_bins.plot ${out}/metabat2_bins.plot
mv ${out}/reassembled_bins.checkm ${out}/work_files

echo " "
echo "########################################################################################################"
echo "########################      BINNING PIPELINE SUCCESSFULLY FINISHED!!!         ########################"
echo "########################################################################################################"
echo " "




