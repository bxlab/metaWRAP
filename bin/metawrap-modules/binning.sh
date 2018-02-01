#!/bin/bash -l

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
	echo "Note: you must chose at least one binning method. Or all at once!"
	echo "Options:"
	echo ""
	echo "	-a STR          metagenomic assembly file"
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo ""
        echo "	--metabat2      bin contigs with metaBAT2"
	echo "	--maxbin2	bin contigs with MaxBin2"
	echo "	--concoct	bin contigs with CONCOCT (warning: this one is slow...)"
	echo "	--run-checkm	immediately run CheckM on the bin results"
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }
run_checkm () {
	comm "Running CheckM on ${1} bins"
	mkdir ${i}.tmp
	checkm lineage_wf -x fa ${1} ${1}.checkm -t $threads --tmpdir ${i}.tmp
	if [[ ! -s ${1}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r ${i}.tmp
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
source config-metawrap

# default params
threads=1; out=false; ASSEMBLY=false
# long options defaults
metabat2=false; maxbin2=false; concoct=false; checkm=false

# load in params
OPTS=`getopt -o ht:o:a: --long help,metabat2,maxbin2,concoct,run-checkm -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
                -o) out=$2; shift 2;;
                -a) ASSEMBLY=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
		--metabat2) metabat2=true; shift 1;;
		--maxbin2) maxbin2=true; shift 1;;
		--concoct) concoct=true; shift 1;;
		--run-checkm) checkm=true; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ $out = false ] || [ $ASSEMBLY = false ] ; then 
	comm "Non-optional parameters -a and/or -o were not entered"
	help_message; exit 1
fi

#check if the assembly file exists
if [ ! -s $ASSEMBLY ]; then error "$ASSEMBLY does not exist. Exiting..."; fi

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


# Make sure at least one binning method was chosen
if [ $metabat2 = false ] && [ $maxbin2 = false ] && [ $concoct = false ]; then
	error "You must select at least one binning method: --metabat2, --maxbin2, --concoct"
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
else echo "Warning: $out already exists."
fi

if [ ! -d ${out}/work_files ]; then mkdir ${out}/work_files; fi

comm "making copy of assembly file $ASSEMBLY"
cp $ASSEMBLY ${out}/work_files/assembly.fa
tmp=${ASSEMBLY##*/}
sample=${tmp%.*}



# Index the assembly
comm "Indexing assembly file"
bwa index ${out}/work_files/assembly.fa
if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi

echo -e "sample\tsample_size\tmean\tstdev" > ${out}/insert_sizes.txt
# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	if [[ $num == *"_1.fastq"* ]]; then 
		reads_1=$num
		reads_2=${num%_*}_2.fastq
		if [ ! -s $reads_1 ]; then error "$reads_1 does not exist. Exiting..."; fi
		if [ ! -s $reads_2 ]; then error "$reads_2 does not exist. Exiting..."; fi

		tmp=${reads_1##*/}
		sample=${tmp%_*}
		
		if [[ ! -f ${out}/work_files/${sample}.bam ]]; then
			comm "Aligning $reads_1 and $reads_2 back to assembly, sorting the alignment, and gathering statistics on insert lengths"
			echo -n -e "${sample}\t" >> ${out}/insert_sizes.txt
			bwa mem -t $threads ${out}/work_files/assembly.fa $reads_1 $reads_2\
			| tee >( awk '{ if ($9 > 0) { N+=1; S+=$9; S2+=$9*$9 }} END { M=S/N; print ""N"\t "M"\t "sqrt ((S2-M*M*N)/(N-1))}'\
			>> ${out}/insert_sizes.txt ) \
			| samtools view -bS - | samtools sort -T ${out}/work_files/tmp-samtools -@ $threads -O bam \
			-o ${out}/work_files/${sample}.bam -
			
			if [[ $? -ne 0 ]]; then error "Something went wrong with aligning/sorting the reads to the assembly!"; fi
		else
			comm "skipping aligning $sample reads to assembly because ${out}/work_files/${sample}.bam already exists."
		fi
	fi
done




if [ $metabat2 = true ]; then
	########################################################################################################
	########################                   RUNNING METABAT2                     ########################
	########################################################################################################
	announcement "RUNNING METABAT2"

	comm "making contig depth file..."	
	jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/metabat2_depth.txt ${out}/work_files/*.bam
	if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi

	comm "Starting binning with metaBAT2..."
	metabat2 -i ${out}/work_files/assembly.fa -a ${out}/work_files/metabat2_depth.txt\
	 -o ${out}/metabat2_bins/bin -m 1500 -t $threads --unbinned
	if [[ $? -ne 0 ]]; then error "Something went wrong with running MetaBAT. Exiting"; fi
	comm "metaBAT2 finished successfully, and found $(ls -l ${out}/metabat2_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/metabat2_bins
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
	
	comm "Starting binning with MaxBin2..."
	mkdir ${out}/work_files/maxbin2_out
	run_MaxBin.pl -contig ${out}/work_files/assembly.fa -markerset 40 -thread $threads\
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
        ########################                   RUNNING CONCOCT                      ########################
        ########################################################################################################
        announcement "RUNNING CONCOCT"

        comm "making contig depth file..."
        jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/tmp --noIntraDepthVariance ${out}/work_files/*.bam
        if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi
        grep -v totalAvgDepth ${out}/work_files/tmp | cut -f 1,4- > ${out}/work_files/concoct_depth.txt
        rm ${out}/work_files/tmp

        comm "Starting binning with CONCOCT. This may take a while..."
        # I have to do some directory changing because CONCOCT dumps all files into current directory...
        home=$(pwd)
        mkdir ${out}/work_files/concoct_out
        cd ${out}/work_files/concoct_out

        concoct --coverage_file ${home}/${out}/work_files/concoct_depth.txt --composition_file ${home}/${out}/work_files/assembly.fa
        if [[ $? -ne 0 ]]; then error "Something went wrong with binning with CONCOCT. Exiting."; fi
        cd $home

        comm "splitting contigs into bins"
        mkdir ${out}/concoct_bins
        ${SOFT}/split_concoct_bins.py ${out}/work_files/concoct_out/clustering_gt1000.csv ${out}/work_files/assembly.fa ${out}/concoct_bins
        comm "CONCOCT finished successfully, and found $(ls -l ${out}/concoct_bins | grep .fa | wc -l) bins!"

	if [ $checkm = true ]; then
		run_checkm ${out}/concoct_bins
	fi
fi

rm ${out}/work_files/*bam

########################################################################################################
########################      BINNING PIPELINE SUCCESSFULLY FINISHED!!!         ########################
########################################################################################################
announcement "BINNING PIPELINE SUCCESSFULLY FINISHED!!!"

