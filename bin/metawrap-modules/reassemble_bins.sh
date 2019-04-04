#!/usr/bin/env bash

##############################################################################################################################################################
#
# Improves a set of bins by aligning reads back to the bins and reassembling them.
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: metaWRAP reassemble_bins [options] -o output_dir -b bin_folder -1 reads_1.fastq -2 reads_2.fastq"
	echo ""
	echo "Options:"
	echo ""
	echo "	-b STR		folder with metagenomic bins"
	echo "	-o STR		output directory"
	echo "	-1 STR          forward reads to use for reassembly"
	echo "	-2 STR          reverse reads to use for reassembly"
	echo ""
	echo "	-t INT		number of threads (default=1)"
	echo "	-m INT		memory in GB (default=40)"
	echo "	-c INT		minimum desired bin completion % (default=70)"
	echo "	-x INT		maximum desired bin contamination % (default=10)"
	echo "	-l INT		minimum contig length to be included in reassembly (default=500)"
	echo ""
	echo "	--strict-cut-off	maximum allowed SNPs for strict read mapping (default=2)"
	echo "	--permissive-cut-off	maximum allowed SNPs for permissive read mapping (default=5)"
	echo "	--skip-checkm		dont run CheckM to assess bins"
	echo "	--parallel		run spades reassembly in parallel, but only using 1 thread per bin"
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }

# these functions are for parallelizing the reassembly
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    "$@" 
    printf '%.3d' $? >&3
    )&
}

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
config_file=$(which config-metawrap)
source $config_file


# default params
threads=1; mem=40; comp=70; cont=10; len=500
bins=None; f_reads=None; r_reads=None; out=None
# long options defaults
strict_max=2; permissive_max=5
run_checkm=true
run_parallel=false
# load in params
OPTS=`getopt -o ht:m:o:x:c:l:b:1:2: --long help,parallel,skip-checkm,strict-cut-off,permissive-cut-off -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
                -m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
		-x) cont=$2; shift 2;;
		-c) comp=$2; shift 2;;
		-b) bins=$2; shift 2;;
		-l) len=$2; shift 2;;
		-1) f_reads=$2; shift 2;;
		-2) r_reads=$2; shift 2;;
                -h | --help) help_message; exit 0; shift 1;;
		--strict-cut-off) strict_max=$2; shift 2;;
		--permissive-cut-off) permissive_max=$2; shift 2;;
		--skip-checkm) run_checkm=false; shift 1;;
		--parallel) run_parallel=true; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ $out = None ] || [  $bins = None ] || [ $f_reads = None ] || [ $r_reads = None ] ; then 
	comm "Some non-optional parameters were not entered"
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi

########################################################################################################
########################               BEGIN REASSEMBLY PIPELINE!               ########################
########################################################################################################
announcement "BEGIN PIPELINE!"
comm "setting up output folder and copying over bins..."
if [ ! -d $out ]; then
        mkdir $out;
else
        warning "Warning: $out already exists!"
fi

if [ -d ${out}/original_bins ]; then rm -r ${out}/original_bins; fi
cp -r $bins ${out}/original_bins
if [ ! -d ${out}/binned_assembly ]; then mkdir ${out}/binned_assembly; fi

# combinde the bins into one big assembly file
if [ -s ${out}/binned_assembly/assembly.fa ]; then rm ${out}/binned_assembly/assembly.fa; fi
for i in $(ls ${out}/original_bins); do cat ${out}/original_bins/$i >> ${out}/binned_assembly/assembly.fa; done


########################################################################################################
########################        RECRUITING READS TO BINS FOR REASSEMBLY         ########################
########################################################################################################
announcement "RECRUITING READS TO BINS FOR REASSEMBLY"

ulimit -n 10000
if [[ $? -ne 0 ]]; then
	ULIMIT=$(ulimit -n)
	warning "Your operating system will allow you to process up to $ULIMIT files at a time. If this is number is less than 4X times of the number of bins you are reassembling, you will likely get an error. Try re-assembling a smaller number of bins."
fi


if [[ ! -s ${out}/binned_assembly/assembly.fa.amb ]]; then
	comm "Indexing the assembly"
	bwa index ${out}/binned_assembly/assembly.fa
	if [[ $? -ne 0 ]]; then error "BWA failed to index $i"; fi

	if [ -d ${out}/reads_for_reassembly ]; then rm -r ${out}/reads_for_reassembly; fi
	mkdir ${out}/reads_for_reassembly

	comm "Aligning all reads back to entire assembly and splitting reads into individual fastq files based on their bin membership"
	bwa mem -t $threads ${out}/binned_assembly/assembly.fa $f_reads $r_reads \
	 | ${SOFT}/filter_reads_for_bin_reassembly.py ${out}/original_bins ${out}/reads_for_reassembly $strict_max $permissive_max
	if [[ $? -ne 0 ]]; then error "Something went wrong with pulling out reads for reassembly..."; fi
else
	comm "WARNING: Looks like the assembly was already indexed. Skipping indexing, and also skipping splitting the reads, because it is assumed you already got to this stage. Will proceed directly to assembly. This is because your output folder $out already has outputs from previous runs. If this is not what you intended, re-run the module with a new output folder, or clear $out."
fi

########################################################################################################
########################             REASSEMBLING BINS WITH SPADES              ########################
########################################################################################################
announcement "REASSEMBLING BINS WITH SPADES"
mkdir ${out}/reassemblies

assemble () {
	bin_name=${1%_*}
	if [[ -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]]; then
		comm "Looks like $bin_name was already re-assembled. Skipping..."
	else
		tmp_dir=${out}/reassemblies/${bin_name}.tmp
		mkdir $tmp_dir
		comm "NOW REASSEMBLING ${bin_name}"
		spades.py -t $2 -m $mem --tmp $tmp_dir --careful \
		--untrusted-contigs ${out}/original_bins/${bin_name%.*}.fa \
		-1 ${out}/reads_for_reassembly/${1%_*}_1.fastq \
		-2 ${out}/reads_for_reassembly/${1%_*}_2.fastq \
		-o ${out}/reassemblies/${bin_name}
		
		if [[ ! -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]]; then
	                warning "Something went wrong with reassembling ${bin_name}"
		else 
			comm "${bin_name} was reassembled successfully!"
			rm -r $tmp_dir
		fi
	fi
}


if [ "$run_parallel" = true ]; then
	open_sem $threads
	for i in $(ls ${out}/reads_for_reassembly/ | grep _1.fastq); do 
		run_with_lock assemble $i 1
	done

	wait
	sleep 1
	comm "all assemblies complete"
fi

if [ "$run_parallel" = false ]; then
	for i in $(ls ${out}/reads_for_reassembly/ | grep _1.fastq); do
		assemble $i $threads
	done

	comm "all assemblies complete"
fi



# removing short contigs and placing reassemblies in the final folder
comm "Finalizing reassemblies"
mkdir ${out}/reassembled_bins
for i in $( ls ${out}/reassemblies/ ); do
	spades_folder=${out}/reassemblies/$i
	bin_name=${spades_folder##*/}
	
	#remove shortest contigs (probably artifacts...)
	if [ -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]; then
		${SOFT}/rm_short_contigs.py $len\
		 ${out}/reassemblies/${bin_name}/scaffolds.fasta\
		 > ${out}/reassemblies/${bin_name}/long_scaffolds.fasta

		if [ -s ${out}/reassemblies/${bin_name}/long_scaffolds.fasta ]; then
			echo "$bin_name was reassembled! Processing..."
			mv ${out}/reassemblies/${bin_name}/long_scaffolds.fasta\
			${out}/reassembled_bins/${bin_name}.fa
		else
			comm "$bin_name was reassembled, but did not yeild contigs $len bp. It is possible there were not enough reads."
		fi
	else
		comm "$bin_name was not successfully reassembled. It is possible there were not enough reads."
	fi
done


if [[ $(ls ${out}/reassembled_bins/ | wc -l) -lt 1 ]]; then
	error "None of the bins were successfully reassembled. ${out}/reassembled_bins/ is empty."
else
	comm "Looks like the reassemblies went well. Now to see if they made the bins better or worse..."
fi



if [ "$run_checkm" = true ]; then
	########################################################################################################
	########################             RUN CHECKM ON REASSEMBLED BINS             ########################
	########################################################################################################
	announcement "RUN CHECKM ON REASSEMBLED BINS"

	# determine --pplacer_threads count. It is either the max thread count or RAM/4, whichever is higher
	ram_max=$(($mem / 40))
	if (( $ram_max < $threads )); then
		p_threads=$ram_max
	else
		p_threads=$threads
	fi
	comm "There is $mem RAM and $threads threads available, and each pplacer thread uses ~40GB, so I will use $p_threads threads for pplacer"

	# copy over original bins
	for base in $( ls ${out}/original_bins/ | grep "\.fa$" ); do 
		i=${out}/original_bins/$base
		cp $i ${out}/reassembled_bins/${base%.*}.orig.fa
	done

	comm "Running CheckM on best bins (reassembled and original)"
	if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
	mkdir ${out}/tmp
	checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
	if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
	rm -r ${out}/tmp


	########################################################################################################
        ########################          FINDING THE BEST VERSION OF EACH BIN          ########################
	########################################################################################################
	announcement "FINDING THE BEST VERSION OF EACH BIN"

	if [ ! -d ${out}/reassembled_best_bins ]; then mkdir ${out}/reassembled_best_bins; fi
	for i in $(${SOFT}/choose_best_bin.py ${out}/reassembled_bins.stats $comp $cont); do 
		echo "Copying best bin: $i"
		cp ${out}/reassembled_bins/${i}.fa ${out}/reassembled_best_bins 
	done
	
	o=$(ls -l ${out}/reassembled_best_bins | grep orig | wc -l)
	s=$(ls -l ${out}/reassembled_best_bins | grep strict | wc -l)
	p=$(ls -l ${out}/reassembled_best_bins | grep permissive | wc -l)

	announcement "Reassembly results are in! $s bins were improved with 'strict' reassembly, $p bins were improved with 'permissive' reassembly, and $o bins were not improved by any reassembly, and thus will stay the same."
	
	if [[ $(ls ${out}/reassembled_best_bins | wc -l) -gt 0 ]]; then 
		comm "Seems that the reassembly went well. You will find the final, best, reassembled bins in ${out}/reassembled_bins, and all intermediate files in ${out}/work_files (which we recomend you delete to save space after you confirm that the pipeline worked)"
		mkdir ${out}/work_files
		mv ${out}/reassembled_bins ${out}/work_files/
		mv ${out}/reassembled_bins.checkm ${out}/work_files/
		mv ${out}/reassembled_bins.stats ${out}/work_files/
		mv ${out}/reads_for_reassembly ${out}/work_files/
		mv ${out}/binned_assembly ${out}/work_files/
		mv ${out}/reassemblies ${out}/work_files/
		#rm -r ${out}/original_bins
		mv ${out}/reassembled_best_bins ${out}/reassembled_bins 
	else
		error "there are no good bins found in ${out}/reassembled_best_bins - something went wrong with choosing the best bins between the reassemblies."
	fi


	comm "Re-running CheckM on the best reasembled bins."
	if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
	mkdir ${out}/tmp
        checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
        if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r ${out}/tmp
        comm "Finalizing CheckM stats..."
        ${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
        if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi

        comm "Making CheckM plot of ${out}/reassembled_bins bins"
        checkm bin_qa_plot -x fa ${out}/reassembled_bins.checkm ${out}/reassembled_bins ${out}/reassembled_bins.plot
        if [[ ! -s ${out}/reassembled_bins.plot/bin_qa_plot.png ]]; then warning "Something went wrong with making the CheckM plot. Exiting."; fi
        mv ${out}/reassembled_bins.plot/bin_qa_plot.png ${out}/reassembled_bins.png
        rm -r ${out}/reassembled_bins.plot
	
	comm "you will find the info on the final reassembled bins in ${out}/reassembled_bins.stats, and a figure summarizing it in ${out}/reassembled_bins.png"

	comm "making reassembly N50, compleiton, and contamination summary plots."
	head -n 1 ${out}/work_files/reassembled_bins.stats > ${out}/original_bins.stats
	grep orig ${out}/work_files/reassembled_bins.stats >> ${out}/original_bins.stats
	${SOFT}/plot_reassembly.py $out $comp $cont ${out}/reassembled_bins.stats ${out}/original_bins.stats
	if [[ $? -ne 0 ]]; then error "Something went wrong with plotting the reassembly summary plots. Exiting..."; fi
fi

comm "you will find the final bins in ${out}/reassembled_bins"
########################################################################################################
########################    REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!        ########################
########################################################################################################
announcement "BIN REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!"

