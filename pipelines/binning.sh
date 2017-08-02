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
	echo ""
	echo "Usage: metaWRAP binning [options] -a assembly.fa -o output_dir readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]"
	echo "Options:"
	echo ""
	echo "	-a STR          metagenomic assembly file"
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT          memory in GB (default=4)"
	echo ""
	echo "	--skip-reassembly	dont reassemble (saves a lot of time)"
	echo "	--skip-checkm		dont run CheckM to assess bins"
	echo "	--checkm-good-bins	re-run CheckM on only good bins (completion>20% and contamination<10%) [Note: pre-reassembly]"
        echo "	--checkm-best-bins      re-run CheckM on the best ressembled version of each bin (good bins only) [Note: post-reassembly]"
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }


# runs CheckM mini-pipeline on a single folder of bins
run_checkm () {
	if [[ -d ${1}.checkm ]]; then rm -r ${1}.checkm; fi
        comm "Running CheckM on $1 bins"
        checkm lineage_wf -x fa $1 ${1}.checkm -t $threads
        if [[ ! -s ${1}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
        comm "Making CheckM plot of $1 bins"
        checkm bin_qa_plot -x fa ${1}.checkm $1 ${1}.plot
        if [[ ! -s ${1}.plot/bin_qa_plot.png ]]; then error "Something went wrong with making the CheckM plot. Exiting."; fi

        comm "Finalizing CheckM stats and plots..."
        ${SOFT}/summarize_checkm.py ${1}.checkm/storage/bin_stats_ext.tsv > ${1}.stats
        mv ${1}.plot/bin_qa_plot.png ${1}.png
        rm -r ${1}.plot
        mv ${1}.checkm ${1%/*}/work_files
}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
source config-metawrap

# default params
threads=1; mem=4; out="false"; ASSEMBLY="false";
# long options defaults
reassemble=true; run_checkm=true; rerun_checkm_on_good_bins=false

# load in params
OPTS=`getopt -o ht:m:o:a: --long help,skip-reassembly,skip-checkm,checkm-good-bins,checkm-best-bins -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
                -m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
                -a) ASSEMBLY=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
		--skip-reassembly) reassemble=false; shift 1;;
		--skip-checkm) run_checkm=false; shift 1;;
		--checkm-good-bins) rerun_checkm_on_good_bins=true; shift 1;;
                --checkm-best-bins) rerun_checkm_on_best_bins=true; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$ASSEMBLY" = "false" ] ; then 
	comm "Non-optional parameters -a and/or -o were not entered"
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


########################################################################################################
########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################
########################################################################################################
announcement "ALIGNING READS TO MAKE COVERAGE FILES"

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

echo -e "sample\tsample_size\tmean\tstdev" > ${out}/insert_sizes.txt
# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	if [[ $num == *"_1.fastq"* ]]; then 
		reads_1=$num
		reads_2=${num%_*}_2.fastq
		if [ "$reads_2" = "false" ]; then error "$reads_2 does not exist. Exiting..."; fi

		tmp=${reads_1##*/}
		sample=${tmp%_*}

		comm "Aligning $reads_1 and $reads_2 back to assembly"
		bwa mem -t $threads ${out}/work_files/alignments/assembly.fa $reads_1 $reads_2 | samtools view -bS - \
		 | samtools sort -T tmp-samtools -@ $threads -O bam -o ${out}/work_files/alignments/${sample}.bam -
		if [[ ! -s ${out}/work_files/alignments/${sample}.bam ]]; then error "Something went wrong with aligning reads to contigs. Exiting."; fi
		
		comm "Saving paired read insert length average and stdev for sample $sample in ${out}/insert_sizes.txt"
		echo -n -e "${sample}\t" >> ${out}/insert_sizes.txt
		samtools view ${out}/work_files/alignments/${sample}.bam | head -n 10000 |\
		 awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print ""N"\t "M"\t "sqrt ((S2-M*M*N)/(N-1))}'\
		 >> ${out}/insert_sizes.txt
	fi
done


########################################################################################################
########################                   RUNNING METABAT2                     ########################
########################################################################################################
announcement "RUNNING METABAT2"

jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/depth.txt ${out}/work_files/alignments/*.bam
metabat2 -i ${out}/work_files/alignments/assembly.fa -a ${out}/work_files/depth.txt\
 -o ${out}/metabat2_bins/bin -m 1500 --unbinned

if [[ ! -s ${out}/metabat2_bins/bin.unbinned.fa ]]; then error "Something went wrong with running MetaBAT. Exiting"; fi


if [ "$run_checkm" = true ]; then
	########################################################################################################
	########################              RUN CHECKM ON ORIGINAL BINS               ########################
	########################################################################################################
	announcement "RUN CHECKM ON ORIGINAL BINS"

	comm "Running CheckM on original metaBAT2 bins"
	run_checkm ${out}/metabat2_bins
	num=$(cat ${out}/metabat2_bins.stats | awk '{if ($2>=20 && $2<=100 && $3>=0 && $3<=10) print $1 }' | wc -l)

	if [ $num -lt 1 ]; then
		comm "There were no 'good' bins detected (>20% completion, <10% contamination)"
 	else
		comm "There are $num 'good' bins found! (>20% completion and <10% contamination) Saving into ${out}/good_bins"
		#find good bins and save them into a seperate folder (completion>20% and contamination<10%)
		mkdir ${out}/good_bins
		for i in $(cat ${out}/metabat2_bins.stats | awk '{if ($2>=20 && $2<=100 && $3>=0 && $3<=10) print $1 }'); do
			cp ${out}/metabat2_bins/${i}.fa ${out}/good_bins
		done
	fi


	# re-run CheckM on only the good bins (makes a more appealing figure...)
	if [ "$rerun_checkm_on_good_bins" = true ]; then
		#check number of good bins
		if [ $num -lt 1 ]; then
			comm "Cant run CheckM on good bins - there were no 'good' bins detected (>20% completion, <10% contamination)"
		else
			comm "Re-running CheckM on only good bins for that nice figure..."
			run_checkm ${out}/good_bins	
		fi
	fi
fi


if [ "$reassemble" = false ]; then
	########################################################################################################
	########################      BINNING PIPELINE FINISHED (NO REASSEMBLY)         ########################
	########################################################################################################
	announcement "BINNING PIPELINE FINISHED (NO REASSEMBLY)"
	exit 1
fi



########################################################################################################
########################        RECRUITING READS TO BINS FOR REASSEMBLY         ########################
########################################################################################################
announcement "RECRUITING READS TO BINS FOR REASSEMBLY"

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




########################################################################################################
########################             REASSEMBLING BINS WITH SPADES              ########################
########################################################################################################
announcement "REASSEMBLING BINS WITH SPADES"

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
	
	if [ -s ${out}/work_files/reassembly/reassemblies/${bin_name}/long_scaffolds.fasta ]; then	
		cp ${out}/work_files/reassembly/reassemblies/${bin_name}/long_scaffolds.fasta\
		 ${out}/reassembled_bins/${bin_name}.fa
	fi
done

if [[ ! -s ${out}/reassembled_bins/bin.unbinned.permissive.fa ]]; then
	error "Something went wrong with processing the reassembled bins and placing them into the reassembled_bins folder. Exiting." 
fi

#REMOVING INTERMEDIATE FILES
comm "Removing intermediate files...."
mv ${out}/work_files/reassembly/reassemblies ${out}/work_files/
rm -r ${out}/work_files/reassembly


if [ "$run_checkm" = true ]; then
	########################################################################################################
	########################             RUN CHECKM ON REASSEMBLED BINS             ########################
	########################################################################################################
	announcement "RUN CHECKM ON REASSEMBLED BINS"

	# copy over original bins
	for i in ${out}/metabat2_bins/*.fa; do base=${i##*/}; cp $i ${out}/reassembled_bins/${base%.*}.orig.fa; done

	comm "Running CheckM on reassembled metaBAT bins"
	run_checkm ${out}/reassembled_bins


	########################################################################################################
        ########################          FINDING THE BEST VERSION OF EACH BIN          ########################
	########################################################################################################
	announcement "FINDING THE BEST VERSION OF EACH BIN"

	cat ${out}/reassembled_bins.stats
	echo ""

	mkdir ${out}/reassembled_best_bins
	for i in $(${SOFT}/choose_best_bin.py ${out}/reassembled_bins.stats); do 
		echo "Copying best bin: $i"
		cp ${out}/reassembled_bins/${i}.fa ${out}/reassembled_best_bins; done
	
	if [ "$rerun_checkm_on_best_bins" = true ]; then
		num=$(ls -l ${out}/reassembled_best_bins | grep .fa | wc -l)
		if [ $num -lt 1 ]; then
			rm -r ${out}/reassembled_best_bins
			warning "There are no 'good' bins priduces. Nothing to run CheckM on..."
		else
			comm "$num good bins found! Re-running CheckM on the best reasembled bins."
			run_checkm ${out}/reassembled_best_bins
		fi
	fi
fi


########################################################################################################
########################      BINNING PIPELINE SUCCESSFULLY FINISHED!!!         ########################
########################################################################################################
announcement "BINNING PIPELINE SUCCESSFULLY FINISHED!!!"

