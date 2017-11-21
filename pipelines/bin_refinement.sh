#!/bin/bash -l

##############################################################################################################################################################
#
# This script is meant to be run on the outputs of binning.sh pipeline to analyze the metagenomic bins and arrive at the best possible putative genomes.
# 
# 1) Use binning_refiner.py to refine bins between 2 or 3 binning folders. Because this script takes in a maximum of 3 bin folders, there are 4 possible
# refinement permutaitons possible (1+2, 2+3, 1+3, and 1+2+3). 
# 2) CheckM is run on all the possible bin folders (1-3 original folders plus 1-4 refinement folders). This gives us completion and contamination estimates
# for every bin. Many bins will be present in varying completions and contaminations in different folders.
# 3) Consolidate all the binning folders by finding the best version of each bin, and saving it into a final "best_bins" folder.
# 4) bwa is used to recruit reads back to the bins and SPAdes ressembles them with. CheckM is again used to find the best version of each bin.
# 5) KRAKEN is used to show the taxonomic distribuion of each bin.
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
# 
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: metaWRAP bin_refinement [options] -o output_dir -A bin_folderA [-B bin_folderB -C bin_folderC] [readsA_1.fastq readsA_2.fastq]"
	echo "Note: the contig names in different bin folders must be consistant (must come from the same assembly)."
	echo "Note 2: you may use any number of F and R read files as long as they end in *_1.fastq and *_2.fastq"
	echo ""
	echo "Options:"
	echo ""
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT          memory in GB (default=4)"
	echo "	-A STR		folder with metagenomic bins"
	echo "	-B STR		another folder with metagenomic bins"
	echo "	-C STR		another folder with metagenomic bins" 
	echo "	-c INT		minimum % completion of bins that is acceptable (default=70)"
	echo "	-x INT		maximum % contamination of bins that is acceptable (default=10)"
	echo ""
	echo "	--skip-refinement	dont use binning_refiner to come up with refined bins based on combinations of binner outputs"
	echo "	--skip-checkm		dont run CheckM to assess bins"
	echo "	--skip-consolidation	choose the best version of each bin from all bin refinement iteration"
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
        comm "Finalizing CheckM stats and plots..."
        ${SOFT}/summarize_checkm.py ${1}.checkm/storage/bin_stats_ext.tsv > ${1}.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
}

# makes checkm plot on a folder of bins if run_checkm has already been run
plot_checkm () {
	comm "Making CheckM plot of $1 bins"
	checkm bin_qa_plot -x fa ${1}.checkm $1 ${1}.plot
	if [[ ! -s ${1}.plot/bin_qa_plot.png ]]; then warning "Something went wrong with making the CheckM plot. Exiting."; fi
	mv ${1}.plot/bin_qa_plot.png ${1}.png
	rm -r ${1}.plot
}


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
source config-metawrap

# default params
threads=1; mem=4; out="false"; comp=70; cont=10; x=10; c=70; 
bins1=None; bins2=None; bins3=None
# long options defaults
run_checkm=true; refine=true; cherry_pick=true

# load in params
OPTS=`getopt -o ht:m:o:x:c:A:B:C: --long help,skip-checkm,skip-refinement,skip-consolidation -- "$@"`
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
		-A) bins1=$2; shift 2;;
		-B) bins2=$2; shift 2;;
		-C) bins3=$2; shift 2;;
                -h | --help) help_message; exit 0; shift 1;;
		--skip-checkm) run_checkm=false; shift 1;;
		--skip-refinement) refine=false; shift 1;;
		--skip-consolidation) cherry_pick=false; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done


########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ $out = false ] || [  bins1 = false ] ; then 
	comm "Non-optional parameters -o and/or -1 were not entered"
	help_message; exit 1
fi

# check for at F and R reads unless --skip-reassembly was specified:
if [ $reassemble = true ]; then
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
	
	comm "$num_of_F_read_files forward and $num_of_R_read_files reverse read files detected, which will be used to bin reassembly at the end of the pipeline."
	if [ ! $num_of_F_read_files == $num_of_R_read_files ]; then error "Number of F and R reads must be the same!"; fi
	if [[ $num_of_F_read_files -lt 1 ]]; then error "You must provide at least one set of reads if you want to reassmeble the reads. Otherwise, use the --skip-reassembly option."; fi
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi



########################################################################################################
########################               BEGIN REFINEMENT PIPELINE!               ########################
########################################################################################################
announcement "PEGIN PIPELINE!"
comm "setting up output folder and copything over bins..."
mkdir $out
n_binnings=0
if [[ -d $bins1 ]]; then 
	cp -r $bins1 ${out}/binsA
	n_binnings=$((n_binnings +1))
else
	error "$bins1 is not a valid directory. Exiting."
fi

#copy over the other bin folders if they are specified
if [[ -d $bins2 ]]; then cp -r $bins2 ${out}/binsB; n_binnings=$((n_binnings +1)); fi
if [[ -d $bins3 ]]; then cp -r $bins3 ${out}/binsC; n_binnings=$((n_binnings +1)); fi

# I have to switch directories here - binning_refiner is pretty glitchy..."
home=$(pwd)
cd $out
if [ $refine = true ]; then
	announcement "BEGIN BIN REFINEMENT"	
	if [[ n_binnings -eq 1 ]]; then
		comm "There is only one bin folder, so no refinement of bins possible. Moving on..."
	elif [[ n_binnings -eq 2 ]]; then
		comm "There are two bin folders, so we can consolidate them into a third, more refined bin set."
		${SOFT}/binning_refiner.py -1 binsA -2 binsB
		mv outputs/Refined binsAB
		rm -r outputs
	elif [[ n_binnings -eq 3 ]]; then
		comm "There are three bin folders, so there 4 ways we can refine the bins (A+B, B+C, A+C, A+B+C). Will try all four!"
		
		comm "pricessing A+B+C"
		${SOFT}/binning_refiner.py -1 binsA -2 binsB -3 binsC
		mv outputs/Refined binsABC; rm -r outputs
		
		comm "pricessing A+B"	
		${SOFT}/binning_refiner.py -1 binsA -2 binsB
		mv outputs/Refined binsAB; rm -r outputs
	
		comm "pricessing B+C"
		${SOFT}/binning_refiner.py -1 binsC -2 binsB
		mv outputs/Refined binsBC; rm -r outputs
		
		comm "pricessing A+C"
		${SOFT}/binning_refiner.py -1 binsA -2 binsC
		mv outputs/Refined binsAC; rm -r outputs
	else
		error "Something is off here - somehow there are not 1, 2, or 3 bin folders ($n_binnings)"
	fi
	comm "Bin refinement finished successfully!"
else
	comm "Skipping bin refinement. Will proceed with the $n_binnings bins specified."
fi
	
comm "fix bin naming to .fa convention for consistancy"
for i in $(ls); do for j in $(ls $i | grep .fasta); do mv ${i}/${j} ${i}/${j%.*}.fa; done; done


########################################################################################################
########################              RUN CHECKM ON ALL BIN SETS                ########################
########################################################################################################
if [ "$run_checkm" = true ]; then
	announcement "RUNNING CHECKM ON ALL SETS OF BINS"

	for bin_set in $(ls); do 
		comm "Running CheckM on $bin_set bins"
		checkm lineage_wf -x fa $bin_set ${bin_set}.checkm -t $threads
		if [[ ! -s ${bin_set}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
		${SOFT}/summarize_checkm.py ${bin_set}.checkm/storage/bin_stats_ext.tsv $bin_set | (read -r; printf "%s\n" "$REPLY"; sort) > ${bin_set}.stats
		if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
		rm -r ${bin_set}.checkm

		num=$(cat ${bin_set}.stats | awk -v c="$comp" -v x="$cont" '{if ($2>=c && $2<=100 && $3>=0 && $3<=x) print $1 }' | wc -l)
		comm "There are $num 'good' bins found in $bin_set! (>${comp}% completion and <${cont}% contamination)"
	done
else
	comm "Skipping CheckM. Warning: bin cnosolidation will not be possible."
fi

########################################################################################################
########################               CONSOLIDATE ALL BIN SETS                 ########################
########################################################################################################
if [ $cherry_pick = true ]; then
	announcement "CONSOLIDATING ALL BIN SETS BY CHOOSING THE BEST VERSION OF EACH BIN"
	if [[ n_binnings -eq 1 ]]; then
	        comm "There is only one original bin folder, so no refinement of bins possible. Moving on..."
		best_bin_set=binsA
	elif [[ n_binnings -eq 2 ]] || [[ n_binnings -eq 3 ]]; then
		comm "There are $n_binnings original bin folders, plus the refined bins."
		rm -r binsM binsM.stats
		cp -r binsA binsM; cp binsA.stats binsM.stats
		for bins in $(ls | grep .stats | grep -v binsM); do
			comm "merging $bins and binsM"
			${SOFT}/consolidate_two_sets_of_bins.py binsM ${bins%.*} binsM.stats $bins binsM1 $comp $cont
			if [[ $? -ne 0 ]]; then error "Something went wrong with merging two sets of bins"; fi
			rm -r binsM binsM.stats
			mv binsM1 binsM; mv binsM1.stats binsM.stats
		done
		best_bin_set=binsM
		comm "All the best bins are stored in $best_bin_set"
	else
		error "Something went wrong with determining the number of bin folders... The number was ${n_binnings}. Exiting."
	fi
	
elif [ $cherry_pick = false ]; then
	comm "Skipping bin consolidation. Will try to pick the best binning folder without mixing bins from different sources."
	if [ $run_checkm = false ]; then 
		comm "cannot decide on best bin set because CheckM was not run. Will assume its binsA (first bin set)"
		best_bin_set=binsA
	elif [ $run_checkm = true ]; then
		max=0
		best_bin_set=none
		for bin_set in $(ls | grep .stats); do
			num=$(cat $bin_set | awk -v c="$comp" -v x="$cont" '{if ($2>=c && $2<=100 && $3>=0 && $3<=x) print $1 }' | wc -l)
			comm "There are $num 'good' bins found in ${bin_set%.*}! (>${comp}% completion and <${cont}% contamination)"
			if [ "$num" -gt "$max" ]; then
				max=$num
				best_bin_set=${bin_set%.*}
			fi
		done
		if [[ ! -d $best_bin_set ]]; then error "Something went wrong with deciding on the best bin set. Exiting."; fi
		comm "looks like the best bin set is $best_bin_set"
	else
		error "something is wrong with the run_checkm option (${run_checkm})"
	fi
else
	error "something is wrong with the cherry_pick option (${cherry_pick})"
fi


comm "moving over temp files into work sub-directory (feel free to delete this)"
mkdir refinement_intermediates
for d in $(ls | grep -v refinement_intermediates); do mv $d refinement_intermediates/; done
cp -r refinement_intermediates/${best_bin_set} best_non_reassembled_bins
cp refinement_intermediates/${best_bin_set}.stats best_non_reassembled_bins.stats

comm "You will find the best non-reassembled versions of the bins in ${out}/best_non_reassembled_bins"

########################################################################################################
########################     BIN_REFINEMENT PIPELINE SUCCESSFULLY FINISHED!!!   ########################
########################################################################################################
announcement "BIN_REFINEMENT PIPELINE FINISHED SUCCESSFULLY! (NO REASSEMBLY)"

