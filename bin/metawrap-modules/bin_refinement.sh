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
	echo "Usage: metaWRAP bin_refinement [options] -o output_dir -A bin_folderA [-B bin_folderB -C bin_folderC]"
	echo "Note: the contig names in different bin folders must be consistant (must come from the same assembly)."
	echo ""
	echo "Options:"
	echo ""
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-A STR		folder with metagenomic bins"
	echo "	-B STR		another folder with metagenomic bins"
	echo "	-C STR		another folder with metagenomic bins" 
	echo "	-c INT		minimum % completion of bins [should be >50%] (default=70)"
	echo "	-x INT		maximum % contamination of bins that is acceptable (default=10)"
	echo ""
	echo "	--skip-refinement	dont use binning_refiner to come up with refined bins based on combinations of binner outputs"
	echo "	--skip-checkm		dont run CheckM to assess bins"
	echo "	--skip-consolidation	choose the best version of each bin from all bin refinement iteration"
	echo "	--keep-ambiguous	for contigs that end up in more than one bin, keep them in all bins (default: keeps them only in the best bin)"
	echo "	--remove-ambiguous	for contigs that end up in more than one bin, remove them in all bins (default: keeps them only in the best bin)"
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
        ${SOFT}/summarize_checkm.py ${1}.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort -rn -k2) > ${1}.stats
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
threads=1; out="false"; comp=70; cont=10; x=10; c=70; 
bins1=None; bins2=None; bins3=None
# long options defaults
run_checkm=true; refine=true; cherry_pick=true; dereplicate=partial

# load in params
OPTS=`getopt -o ht:o:x:c:A:B:C: --long help,skip-checkm,skip-refinement,skip-consolidation,keep-ambiguous,remove-ambiguous -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
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
		--keep-ambiguous) dereplicate=false; shift 1;;
		--remove-ambiguous) dereplicate=complete; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [[ $out == false ]] || [[  $bins1 == false ]] ; then 
	comm "Non-optional parameters -o and/or -A were not entered"
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
fi


########################################################################################################
########################               BEGIN REFINEMENT PIPELINE!               ########################
########################################################################################################
announcement "BEGIN PIPELINE!"
comm "setting up output folder and copything over bins..."
if [ ! -d $out ]; then
        mkdir $out
	if [ ! -d $out ]; then error "cannot make $out"; fi
else
        echo "Warning: $out already exists. Cleaning..."
	rm -r ${out}/*
fi

n_binnings=0
if [[ -d $bins1 ]]; then 
	cp -r $bins1 ${out}/binsA
	n_binnings=$((n_binnings +1))
else
	error "$bins1 is not a valid directory. Exiting."
fi

comm "there are $(ls ${out}/binsA | wc -l) bins in binsA"
if [[ $(ls ${out}/binsA | wc -l) -eq 0 ]]; then error "Please provide valid input. Exiting..."; fi
#copy over the other bin folders if they are specified
if [[ -d $bins2 ]]; then 
	cp -r $bins2 ${out}/binsB; n_binnings=$((n_binnings +1))
	comm "there are $(ls ${out}/binsB | wc -l) bins in binsB"
	if [[ $(ls ${out}/binsB | wc -l) -eq 0 ]]; then error "Please provide valid input. Exiting..."; fi
fi
if [[ -d $bins3 ]]; then 
	cp -r $bins3 ${out}/binsC; n_binnings=$((n_binnings +1))
	comm "there are $(ls ${out}/binsC | wc -l) bins in binsC"
	if [[ $(ls ${out}/binsC | wc -l) -eq 0 ]]; then error "Please provide valid input. Exiting..."; fi
fi

comm "There are $n_binnings bin sets!"

# I have to switch directories here - Bin_refiner dumps everything into the current dir"
home=$(pwd)
cd $out

if [ "$refine" == "true" ]; then
	announcement "BEGIN BIN REFINEMENT"	
	if [[ $n_binnings -eq 1 ]]; then
		comm "There is only one bin folder, so no refinement of bins possible. Moving on..."
	elif [[ $n_binnings -eq 2 ]]; then
		comm "There are two bin folders, so we can consolidate them into a third, more refined bin set."
		${SOFT}/binning_refiner.py -1 binsA -2 binsB -o Refined_AB
		mv Refined_AB/Refined binsAB
		if [[ $? -ne 0 ]]; then error "Bin_refiner did not finish correctly. Exiting..."; fi
		rm -r Refined_AB
	elif [[ $n_binnings -eq 3 ]]; then
		comm "There are three bin folders, so there 4 ways we can refine the bins (A+B, B+C, A+C, A+B+C). Will try all four in parallel!"
		
		${SOFT}/binning_refiner.py -1 binsA -2 binsB -3 binsC -o Refined_ABC &
		${SOFT}/binning_refiner.py -1 binsA -2 binsB -o Refined_AB &
		${SOFT}/binning_refiner.py -1 binsC -2 binsB -o Refined_BC &
		${SOFT}/binning_refiner.py -1 binsA -2 binsC -o Refined_AC &
		
		wait
	
		comm "there are $(ls Refined_AB/Refined | wc -l) refined bins in binsAB"
		comm "there are $(ls Refined_BC/Refined | wc -l) refined bins in binsBC"
		comm "there are $(ls Refined_AC/Refined | wc -l) refined bins in binsAC"
		comm "there are $(ls Refined_ABC/Refined | wc -l) refined bins in binsABC"


		mv Refined_ABC/Refined binsABC
		if [[ $? -ne 0 ]]; then error "Bin_refiner did not finish correctly with A+B+C. Exiting..."; fi
		rm -r Refined_ABC
		
		mv Refined_AB/Refined binsAB
		if [[ $? -ne 0 ]]; then error "Bin_refiner did not finish correctly with A+B. Exiting..."; fi
		rm -r Refined_AB
	
		mv Refined_BC/Refined binsBC
		if [[ $? -ne 0 ]]; then error "Bin_refiner did not finish correctly with B+C. Exiting..."; fi
		rm -r Refined_BC
		
		mv Refined_AC/Refined binsAC
		if [[ $? -ne 0 ]]; then error "Bin_refiner did not finish correctly with A+C. Exiting..."; fi
		rm -r Refined_AC
	else
		error "Something is off here - somehow there are not 1, 2, or 3 bin folders ($n_binnings)"
	fi
	comm "Bin refinement finished successfully!"
else
	comm "Skipping bin refinement. Will proceed with the $n_binnings bins specified."
fi
	
comm "fixing bin naming to .fa convention for consistancy..."
for i in $(ls); do for j in $(ls $i | grep .fasta); do mv ${i}/${j} ${i}/${j%.*}.fa; done; done


########################################################################################################
########################              RUN CHECKM ON ALL BIN SETS                ########################
########################################################################################################
if [ "$run_checkm" == "true" ]; then
	announcement "RUNNING CHECKM ON ALL SETS OF BINS"
	for bin_set in $(ls | grep -v tmp); do 
		comm "Running CheckM on $bin_set bins"
		mkdir ${bin_set}.tmp
		checkm lineage_wf -x fa $bin_set ${bin_set}.checkm -t $threads --tmpdir ${bin_set}.tmp
		if [[ ! -s ${bin_set}.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
		${SOFT}/summarize_checkm.py ${bin_set}.checkm/storage/bin_stats_ext.tsv $bin_set | (read -r; printf "%s\n" "$REPLY"; sort) > ${bin_set}.stats
		if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
		rm -r ${bin_set}.checkm; rm -r ${bin_set}.tmp

		num=$(cat ${bin_set}.stats | awk -v c="$comp" -v x="$cont" '{if ($2>=c && $2<=100 && $3>=0 && $3<=x) print $1 }' | wc -l)
		comm "There are $num 'good' bins found in $bin_set! (>${comp}% completion and <${cont}% contamination)"
	done
else
	comm "Skipping CheckM. Warning: bin cnosolidation will not be possible."
fi

########################################################################################################
########################               CONSOLIDATE ALL BIN SETS                 ########################
########################################################################################################
if [ "$cherry_pick" == "true" ]; then
	announcement "CONSOLIDATING ALL BIN SETS BY CHOOSING THE BEST VERSION OF EACH BIN"
	if [[ $n_binnings -eq 1 ]]; then
	        comm "There is only one original bin folder, so no refinement of bins possible. Moving on..."
		best_bin_set=binsA
	elif [[ $n_binnings -eq 2 ]] || [[ $n_binnings -eq 3 ]]; then
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

		if [[ $dereplicate == false ]]; then
			comm "Skipping dereplication of contigs between bins..."
			mv binsM binsO
			mv binsM.stats binsO.stats
		elif [[ $dereplicate == partial ]]; then
			comm "Scanning to find duplicate contigs between bins and only keep them in the best bin..."
			${SOFT}/dereplicate_contigs_in_bins.py binsM.stats binsM binsO
		elif [[ $dereplicate == complete ]]; then
			comm "Scanning to find duplicate contigs between bins and deleting them in all bins..."
			${SOFT}/dereplicate_contigs_in_bins.py binsM.stats binsM binsO remove
		else
			error "there was an error in deciding how to dereplicate contigs"
		fi

		best_bin_set=binsO
	else
		error "Something went wrong with determining the number of bin folders... The number was ${n_binnings}. Exiting."
	fi
	
elif [ "$cherry_pick" == "false" ]; then
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

comm "You will find the best non-reassembled versions of the bins in $best_bin_set"


########################################################################################################
########################               FINALIZING THE REFINED BINS              ########################
########################################################################################################
announcement "FINALIZING THE REFINED BINS"


if [ "$run_checkm" == "true" ] && [ $dereplicate != "false" ]; then
	comm "Re-running CheckM on binsO bins"
	mkdir binsO.tmp
	checkm lineage_wf -x fa binsO binsO.checkm -t $threads --tmpdir binsO.tmp
	if [[ ! -s binsO.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r binsO.tmp
	${SOFT}/summarize_checkm.py binsO.checkm/storage/bin_stats_ext.tsv manual binsM.stats | (read -r; printf "%s\n" "$REPLY"; sort -rn -k2) > binsO.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
	rm -r binsO.checkm
	num=$(cat binsO.stats | awk -v c="$comp" -v x="$cont" '{if ($2>=c && $2<=100 && $3>=0 && $3<=x) print $1 }' | wc -l)
	comm "There are $num 'good' bins found in binsO.checkm! (>${comp}% completion and <${cont}% contamination)"
	
	comm "Removing bins that are inadequate quality..."
	for bin_name in $(cat binsO.stats | grep -v compl | awk -v c="$comp" -v x="$cont" '{if ($2<c || $2>100 || $3<0 || $3>x) print $1 }' | cut -f1); do
		echo "${bin_name} will be removed because it fell below the quality threshhold after de-replication of contigs..."
		rm binsO/${bin_name}.fa
	done
	head -n 1 binsO.stats > binsO.stats.tmp
	cat binsO.stats | awk -v c="$comp" -v x="$cont" '$2>=c && $2<=100 && $3>=0 && $3<=x' >> binsO.stats.tmp
	mv binsO.stats.tmp binsO.stats
	n=$(cat binsO.stats | grep -v comp | wc -l)
	comm "Re-evaluating bin quality after contig de-replication is complete! There are still $n high quality bins."
fi


if [ "$run_checkm" == "true" ]; then
	comm "making completion and contamination ranking plots for all refinement iterations"
	${SOFT}/plot_binning_results.py $(ls | grep ".stats")
	mkdir figures
	mv binning_results.eps figures/intermediate_binning_results.eps
	mv binning_results.png figures/intermediate_binning_results.png
fi


comm "Finalizing bin set naming, and moving over intermediate files into 'work_files'"
if [[ $n_binnings -ne 1 ]]; then
	mkdir work_files
	for f in $(ls | grep -v results | grep -v figures | grep -v work); do
		mv $f work_files/
	done
fi

if [ "${bins1:$((${#bins1}-1)):1}" = "/" ]; then bins1=${bins1%/*}; fi
if [ "${bins2:$((${#bins2}-1)):1}" = "/" ]; then bins2=${bins2%/*}; fi
if [ "${bins3:$((${#bins3}-1)):1}" = "/" ]; then bins3=${bins3%/*}; fi

if [[ $n_binnings -eq 2 ]]; then
	cp -r work_files/binsA ${bins1##*/}
	cp -r work_files/binsB ${bins2##*/}
	cp -r work_files/binsO metaWRAP_bins

	cp work_files/binsA.stats ${bins1##*/}.stats	
	cp work_files/binsB.stats ${bins2##*/}.stats
	cp work_files/binsO.stats metaWRAP_bins.stats
elif [[ $n_binnings -eq 3 ]]; then
	cp -r work_files/binsA ${bins1##*/}
	cp -r work_files/binsB ${bins2##*/}
	cp -r work_files/binsC ${bins3##*/}
	cp -r work_files/binsO metaWRAP_bins

	cp work_files/binsA.stats ${bins1##*/}.stats
	cp work_files/binsB.stats ${bins2##*/}.stats
	cp work_files/binsC.stats ${bins3##*/}.stats
	cp work_files/binsO.stats metaWRAP_bins.stats
else
	error "something went wrong with determining the number of input bin sets when finalizing output names..."
fi

if [ "$run_checkm" == "true" ]; then
        comm "making completion and contamination ranking plots of final outputs"
        ${SOFT}/plot_binning_results.py $comp $cont $(ls | grep ".stats")
	mv binning_results.eps figures/binning_results.eps
	mv binning_results.png figures/binning_results.png
fi

cd $home

########################################################################################################
########################     BIN_REFINEMENT PIPELINE SUCCESSFULLY FINISHED!!!   ########################
########################################################################################################
announcement "BIN_REFINEMENT PIPELINE FINISHED SUCCESSFULLY!"

