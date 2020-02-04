#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script is meant to be run on the outputs of assembly.sh pipeline to split the assembly contigs into metagenomic bins.
# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.
# The more samples, the better the binning. 
#
# The script uses metaBAT2 and/or CONCOCT and/or MaxBin2 to bin the contigs. MetaBAT2 is the defualt due to its speed and great performance,
# but all these binners have their advantages and disadvantages, so it recomended to run the bin_refinement module to QC the bins, get the 
# best bins of all of each method, and to reassembly and refine the final bins. 
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: metaWRAP binning [options] -a assembly.fa -o output_dir readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]"
	echo "Note1: Make sure to provide all your separately replicate read files, not the joined file."
	echo "Note2: You may provide single end or interleaved reads as well with the use of the correct option"
	echo "Note3: If the output already has the .bam alignments files from previous runs, the module will skip re-aligning the reads"
	echo ""
	echo "Options:"
	echo ""
	echo "	-a STR          metagenomic assembly file"
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT		amount of RAM available (default=4)"
	echo "	-l INT		minimum contig length to bin (default=1000bp). Note: metaBAT will default to 1500bp minimum"
	echo ""
        echo "	--metabat2      bin contigs with metaBAT2"
	echo "	--metabat1	bin contigs with the original metaBAT"
	echo "	--maxbin2	bin contigs with MaxBin2"
	echo "	--concoct	bin contigs with CONCOCT"
	echo ""
	echo "	--universal	use universal marker genes instead of bacterial markers in MaxBin2 (improves Archaea binning)"
	echo "	--run-checkm	immediately run CheckM on the bin results (requires 40GB+ of memory)"
	echo "	--single-end	non-paired reads mode (provide *.fastq files)"
	echo "	--interleaved	the input read files contain interleaved paired-end reads"
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }
run_checkm () {
	comm "Running CheckM on ${1} bins"

	# determine --pplacer_threads count. It is either the max thread count or RAM/4, whichever is higher
	ram_max=$(($mem / 40))
	if (( $ram_max < $threads )); then
		p_threads=$ram_max
	else
		p_threads=$threads
	fi
	comm "There is $mem RAM and $threads threads available, and each pplacer thread uses <40GB, so I will use $p_threads threads for pplacer"

	mkdir ${1}.tmp
	checkm lineage_wf -x fa ${1} ${1}.checkm -t $threads --tmpdir ${1}.tmp --pplacer_threads $p_threads
	if [[ ! -s ${1}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r ${1}.tmp
	${SOFT}/summarize_checkm.py ${1}.checkm/storage/bin_stats_ext.tsv ${1##*/}\
	| (read -r; printf "%s\n" "$REPLY"; sort -rn -k2) > ${1}.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
	num=$(cat ${1}.stats | awk '{if ($2>=70 && $2<=100 && $3>=0 && $3<=10) print $1 }' | wc -l)
	comm "There are $num 'good' bins found in ${1}! (>70% completion and <10% contamination)"
}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
config_file=$(which config-metawrap)
source $config_file

# default params
threads=1; mem=4; len=1000; out=false; ASSEMBLY=false
# long options defaults
metabat1=false; metabat2=false; maxbin2=false; concoct=false
checkm=false; read_type=paired
markers=107

# load in params
OPTS=`getopt -o ht:m:o:a:l: --long help,metabat1,metabat2,maxbin2,concoct,run-checkm,single-end,universal,interleaved -- "$@"`
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
		-m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
                -a) ASSEMBLY=$2; shift 2;;
		-l) len=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
		--metabat2) metabat2=true; shift 1;;
		--metabat1) metabat1=true; shift 1;;
		--maxbin2) maxbin2=true; shift 1;;
		--concoct) concoct=true; shift 1;;
		--run-checkm) checkm=true; shift 1;;
		--single-end) read_type=single; shift 1;;
		--interleaved) read_type=interleaved; shift 1;;
		--universal) markers=40; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################
# Make sure at least one binning method was chosen
if [ $metabat2 = false ] && [ $metabat1 = false ] &&[ $maxbin2 = false ] && [ $concoct = false ]; then
	help_message
	error "You must select at least one binning method: --metabat2, --metabat1, --maxbin2, --concoct"
fi

# check if all parameters are entered
if [ $out = false ] || [ $ASSEMBLY = false ] ; then 
	comm "Non-optional parameters -a and/or -o were not entered"
	help_message; exit 1
fi

#check if the assembly file exists
if [ ! -s $ASSEMBLY ]; then error "$ASSEMBLY does not exist. Exiting..."; fi

comm "Entered read type: $read_type"

if [ $read_type = paired ]; then
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
else
	# check for at least one fastq read
	F="no"
	for num in "$@"; do
		if [[ $num == *".fastq" ]]; then F="yes"; fi
	done
	if [ $F = "no" ]; then
		comm "Unable to find read files in format *.fastq (for single-end or interleaved reads)"
		help_message; exit 1
	fi
fi

if [ $read_type = paired ]; then
	#determine number of fastq read files provided:
	num_of_F_read_files=$(for I in "$@"; do echo $I | grep _1.fastq; done | wc -l)
	num_of_R_read_files=$(for I in "$@"; do echo $I | grep _2.fastq; done | wc -l)

	comm "$num_of_F_read_files forward and $num_of_R_read_files reverse read files detected"
	if [ ! $num_of_F_read_files == $num_of_R_read_files ]; then error "Number of F and R reads must be the same!"; fi
fi

if [ $len -lt 1500 ]; then
	metabat_len=1500
else
	metabat_len=$len
fi


# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################


########################################################################################################
########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################
########################################################################################################
announcement "ALIGNING READS TO MAKE COVERAGE FILES"

# setting up the output folder
if [ ! -d $out ]; then mkdir $out;
else 
	echo "Warning: $out already exists."
	rm -r ${out}/*checkm
fi

if [ ! -d ${out}/work_files ]; then mkdir ${out}/work_files; fi

if [ -f ${out}/work_files/assembly.fa ]; then
	comm "Looks like the assembly file is already coppied, but will re-transfer just in case to avoid truncation problems."
	cp $ASSEMBLY ${out}/work_files/assembly.fa
else
	comm "making copy of assembly file $ASSEMBLY"
	cp $ASSEMBLY ${out}/work_files/assembly.fa
fi

tmp=${ASSEMBLY##*/}
sample=${tmp%.*}

# Index the assembly
if [ -f ${out}/work_files/assembly.fa.bwt ]; then
	comm "Looks like there is a index of the assembly already. Skipping..."
else
	comm "Indexing assembly file"
	bwa index ${out}/work_files/assembly.fa
	if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi
fi

# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	# paired end reads
	if [ $read_type = paired ]; then
		if [[ $num == *"_1.fastq"* ]]; then 
			reads_1=$num
			reads_2=${num%_*}_2.fastq
			if [ ! -s $reads_1 ]; then error "$reads_1 does not exist. Exiting..."; fi
			if [ ! -s $reads_2 ]; then error "$reads_2 does not exist. Exiting..."; fi
	
			tmp=${reads_1##*/}
			sample=${tmp%_*}
			
			if [[ ! -f ${out}/work_files/${sample}.bam ]]; then
				comm "Aligning $reads_1 and $reads_2 back to assembly"
				bwa mem -v 1 -t $threads ${out}/work_files/assembly.fa $reads_1 $reads_2 > ${out}/work_files/${sample}.sam
				if [[ $? -ne 0 ]]; then error "Something went wrong with aligning $reads_1 and $reads_2 reads to the assembly. Exiting"; fi

				comm "Sorting the $sample alignment file"
				samtools sort -T ${out}/work_files/tmp-samtools -@ $threads -O BAM -o ${out}/work_files/${sample}.bam ${out}/work_files/${sample}.sam
				if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiging..."; fi
				rm ${out}/work_files/${sample}.sam
			else
				comm "skipping aligning $sample reads to assembly because ${out}/work_files/${sample}.bam already exists."
			fi
		fi

	# single end or interleaved reads
	else
		if [[ $num == *".fastq"* ]]; then
			reads=$num
			if [ ! -s $reads ]; then error "$reads does not exist. Exiting..."; fi
			tmp=${reads##*/}
			sample=${tmp%.*}
			if [[ ! -f ${out}/work_files/${sample}.bam ]]; then
				comm "Aligning $reads back to assembly, and sorting the alignment"
				if [ $read_type = single ]; then
					bwa mem -t $threads ${out}/work_files/assembly.fa $reads > ${out}/work_files/${sample}.sam
					if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
				elif [ $read_type = interleaved ]; then
					bwa mem -v 1 -p -t $threads ${out}/work_files/assembly.fa $reads > ${out}/work_files/${sample}.sam
					if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
				else
					error "something from with the read_type (=$read_type)"
				fi
				
				comm "Sorting the $sample alignment file"
				samtools sort -T ${out}/work_files/tmp-samtools -@ $threads -O BAM -o ${out}/work_files/${sample}.bam ${out}/work_files/${sample}.sam
				if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiging..."; fi
				rm ${out}/work_files/${sample}.sam
			else
				comm "skipping aligning $sample reads to assembly because ${out}/work_files/${sample}.bam already exists."
			fi
		fi
	fi
done


if [ $metabat2 = true ]; then
	########################################################################################################
	########################                   RUNNING METABAT2                     ########################
	########################################################################################################
	announcement "RUNNING METABAT2"

	comm "making contig depth file..."	
	jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/metabat_depth.txt ${out}/work_files/*.bam
	if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi

	comm "Starting binning with metaBAT2..."
	metabat2 -i ${out}/work_files/assembly.fa -a ${out}/work_files/metabat_depth.txt\
	 -o ${out}/metabat2_bins/bin -m $metabat_len -t $threads --unbinned
	if [[ $? -ne 0 ]]; then error "Something went wrong with running MetaBAT2. Exiting"; fi
	comm "metaBAT2 finished successfully, and found $(ls -l ${out}/metabat2_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/metabat2_bins
	fi
fi


if [ $metabat1 = true ]; then
        ########################################################################################################
        ########################                   RUNNING METABAT1                     ########################
        ########################################################################################################
        announcement "RUNNING METABAT1"

        comm "making contig depth file..."
        jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/metabat_depth.txt ${out}/work_files/*.bam
        if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi

        comm "Starting binning with metaBAT1..."
        metabat1 -i ${out}/work_files/assembly.fa -a ${out}/work_files/metabat_depth.txt\
         -o ${out}/metabat1_bins/bin -m $metabat_len --minContigByCorr $len -t $threads --unbinned --superspecific
        if [[ $? -ne 0 ]]; then error "Something went wrong with running MetaBAT1. Exiting"; fi
        comm "metaBAT1 finished successfully, and found $(ls -l ${out}/metabat1_bins | grep .fa | wc -l) bins!"

        if [ $checkm = true ]; then
                run_checkm ${out}/metabat1_bins
        fi
fi



if [ $maxbin2 = true ]; then
        ########################################################################################################
        ########################                   RUNNING MAXBIN2                     ########################
        ########################################################################################################
        announcement "RUNNING MAXBIN2"

	comm "making contig depth file..."
        jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/mb2_master_depth.txt --noIntraDepthVariance ${out}/work_files/*.bam
        if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi

	#calculate total numper of columns
	A=($(head -n 1 ${out}/work_files/mb2_master_depth.txt)) 
	N=${#A[*]}
	
	# split the contig depth file into multiple files
	comm "split master contig depth file into individual files for maxbin2 input"
	if [ -f ${out}/work_files/mb2_abund_list.txt ]; then rm ${out}/work_files/mb2_abund_list.txt; fi
	for i in $(seq 4 $N); do 
		sample=$(head -n 1 ${out}/work_files/mb2_master_depth.txt | cut -f $i)
		echo "processing $sample depth file..."
		grep -v totalAvgDepth ${out}/work_files/mb2_master_depth.txt | cut -f 1,$i > ${out}/work_files/mb2_${sample%.*}.txt
		if [[ $out == /* ]]; then
			echo ${out}/work_files/mb2_${sample%.*}.txt >> ${out}/work_files/mb2_abund_list.txt
		else
			echo $(pwd)/${out}/work_files/mb2_${sample%.*}.txt >> ${out}/work_files/mb2_abund_list.txt
		fi
	done

	run_MaxBin.pl
	if [[ $? -ne 0 ]]; then
		comm "looks like our default perl libraries are not the conda ones. Manually setting perl5 library directory"
        	conda_path=$(which metawrap)
		echo "metawrap path: $conda_path"
		conda_path=${conda_path%/*}
		if [ $(echo -n $conda_path | tail -c 1) = "/" ]; then conda_path=${conda_path%/*}; fi
		conda_path=${conda_path%/*}
		if [ ! -d ${conda_path}/lib/perl5/site_perl/5.22.0 ]; then 
			error "${conda_path}/lib/perl5/site_perl/5.22.0 does not exixt. Cannot set manual path to perl5 libraries. Exiting..."
		fi
	
	        perl_libs=${conda_path}/lib/perl5/site_perl/5.22.0
	        echo "Will use perl5 libraries located in $perl_libs - hopefully they are there. Install Perl in the conda environment if the directory does not exist."
		export PERL5LIB="$perl_libs"
	fi

	
	comm "Starting binning with MaxBin2..."
	mkdir ${out}/work_files/maxbin2_out
	run_MaxBin.pl -contig ${out}/work_files/assembly.fa -markerset $markers -thread $threads -min_contig_length $len\
	-out ${out}/work_files/maxbin2_out/bin \
	-abund_list ${out}/work_files/mb2_abund_list.txt

	if [[ $? -ne 0 ]]; then error "Something went wrong with running MaxBin2. Exiting."; fi
	if [[ $(ls ${out}/work_files/maxbin2_out/ | grep bin | grep .fasta | wc -l) -lt 1 ]]; then error "MaxBin2 did not pruduce a single bin. Something went wrong. Exiting."; fi

	mkdir ${out}/maxbin2_bins
	N=0
	for i in $(ls ${out}/work_files/maxbin2_out/ | grep bin | grep .fasta); do
		cp ${out}/work_files/maxbin2_out/$i ${out}/maxbin2_bins/bin.${N}.fa
		N=$((N + 1))
	done
	comm "MaxBin2 finished successfully, and found $(ls -l ${out}/maxbin2_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/maxbin2_bins
	fi
fi

if [ $concoct = true ]; then
	########################################################################################################
	########################                    RUNNING CONCOCT                     ########################
	########################################################################################################
        announcement "RUNNING CONCOCT"

	if [[ ! -s ${out}/work_files/concoct_depth.txt ]]; then
		comm "indexing .bam alignment files..."
		for FILE in ${out}/work_files/*.bam; do
			echo $FILE
			samtools index -@ $threads -b $FILE
		done

        	comm "cutting up contigs into 10kb fragments for CONCOCT..."
		cut_up_fasta.py ${out}/work_files/assembly.fa -c 10000 --merge_last -b ${out}/work_files/assembly_10K.bed -o 0 > ${out}/work_files/assembly_10K.fa
        	if [[ $? -ne 0 ]]; then error "Something went wrong with cutting up contigs. Exiting."; fi

		comm "estimating contig fragment coverage..."	
		CMD="concoct_coverage_table.py ${out}/work_files/assembly_10K.bed ${out}/work_files/*.bam > ${out}/work_files/concoct_depth.txt"
		$(eval $CMD)
		if [[ $? -ne 0 ]]; then error "Something went wrong with estimating fragment abundance. Exiting..."; fi
	else
		comm "looks like contig coverage was already estimated... skipping"
	fi


	concoct -h
        if [[ $? -ne 0 ]]; then
                comm "Looks like our environment has a faulty libgslcblas.so link. Will try to manually create symlink in conda environment"
                conda_path=$(which concoct)
                conda_path=${conda_path%/*}
                if [ $(echo -n $conda_path | tail -c 1) = "/" ]; then conda_path=${conda_path%/*}; fi
                conda_path=${conda_path%/*}
		echo "conda path: $conda_path"
                if [ ! -s ${conda_path}/lib/libgslcblas.so ]; then
                        error "${conda_path}/lib/libgslcblas.so does not exixt. Cannot set libgslcblas.so.0 symlink. Please make sure that concoct is working before re-trying. Exiting..."
                fi

		echo "Creating symlink ${conda_path}/lib/libgslcblas.so.0 from ${conda_path}/lib/libgslcblas.so"
		rm ${conda_path}/lib/libgslcblas.so.0
		ln ${conda_path}/lib/libgslcblas.so ${conda_path}/lib/libgslcblas.so.0
	fi


        comm "Starting binning with CONCOCT..."
        mkdir ${out}/work_files/concoct_out

	concoct -l $len -t $threads \
		--coverage_file ${out}/work_files/concoct_depth.txt \
		--composition_file ${out}/work_files/assembly_10K.fa \
		-b ${out}/work_files/concoct_out

	if [[ $? -ne 0 ]]; then error "Something went wrong with binning with CONCOCT. Exiting..."; fi

	comm "merging 10kb fragments back into contigs"
	merge_cutup_clustering.py ${out}/work_files/concoct_out/clustering_gt${len}.csv > ${out}/work_files/concoct_out/clustering_gt${len}_merged.csv
	if [[ $? -ne 0 ]]; then error "Something went wrong with merging fragments. Exiting..."; fi

        comm "splitting contigs into bins"
        mkdir ${out}/concoct_bins
        ${SOFT}/split_concoct_bins.py ${out}/work_files/concoct_out/clustering_gt${len}_merged.csv ${out}/work_files/assembly.fa ${out}/concoct_bins
        comm "CONCOCT finished successfully, and found $(ls -l ${out}/concoct_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/concoct_bins
	fi
fi

#comm "cleaning up *.bam to save space..."
#rm ${out}/work_files/*bam

########################################################################################################
########################      BINNING PIPELINE SUCCESSFULLY FINISHED!!!         ########################
########################################################################################################
announcement "BINNING PIPELINE SUCCESSFULLY FINISHED!!!"

