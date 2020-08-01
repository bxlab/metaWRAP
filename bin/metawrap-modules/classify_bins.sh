#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script classifies all contigs in a set of bins by aligning them to the NCBI database with MEGABLAST, pruning the resulting hits, and assigning final
# taxonomy with taxator-k. The consensus taxonomy of each bin is called by contructing aweighted consensus tree, and traversing the tree by the most
# likely path.
#
# For questions, bugs, and suggestions, contact German Uritskiy at guritsk1@jhu.edu.
# 
##############################################################################################################################################################



help_message () {
	echo ""
	echo "Usage: metaWRAP classify_bins [options] -b bin_folder -o output_dir"
	echo "Options:"
	echo ""
	echo "	-b STR		folder with the bins to be classified (in fasta format)"
	echo "	-o STR		output directory"
	echo "	-t INT          number of threads"
	echo ""
	echo "";}

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }


########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################


# setting scripts and databases from config file (should be in same folder as main script)
config_file=$(which config-metawrap)
source $config_file


# Set defaults
threads=1; out="false"; bin_folder="false"

# load in params
OPTS=`getopt -o ht:o:b: --long help -- "$@"`
# make sure the params are entered correctly
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
                -o) out=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
                -b) bin_folder=$2; shift 2;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done



######################################################################################################o#
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

# check if all parameters are entered
if [ "$out" = "false" ] || [ "$bin_folder" = "false" ]; then 
	help_message; exit 1
fi

# Checks for correctly configures meta-scripts folder
if [ ! -s $SOFT/sort_contigs.py ]; then
	error "The folder $SOFT doesnt exist. Please make sure the meta-scripts folder path is correctly set in config.sh file"
fi

#  Checks for NCBI_nt database for BLAST
if [ ! -f "${BLASTDB}/nt.00.nhd" ]; then
	error "The file ${BLASTDB}/nt.00.nhd doesnt exist, which likely means that you havent set the correct path to your NCBI_nt database or\
	 havent downloaded it. Please look and the the script and manually set the correct path to NCBI_nt, and follow these steps to download the\
	  database: wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz; for a in nt.*.tar.gz; do tar xzf $a; done"
fi

#  Checks for NCBI_tax database for BLAST
if [ ! -f "${TAXDUMP}/names.dmp" ]; then
	error "The file ${TAXDUMP}/citations.dmp doesnt exist, which likely means that you havent set the correct path to your NCBI_tax database,\
	 or havent downloaded it yet. Please look and the the script and manually set the correct path to NCBI_tax, and follow these steps to download \
	 the database: wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz; tar -xvf taxdump.tar.gz;"
fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################



########################################################################################################
########################       ALIGN CONTIGS TO DATABASE WITH MEGABLAST         ########################
########################################################################################################
announcement "ALIGN CONTIGS TO DATABASE WITH MEGABLAST"

comm "setting up ouput folder $out and merging contigs from all bins..."
if [ ! -d $out ]; then
        mkdir $out;
else
        echo "Warning: $out already exists."
fi

if [[ -s ${out}/all_contigs.fa ]]; then rm ${out}/all_contigs.fa; fi
for f in $(ls $bin_folder); do cat ${bin_folder}/${f} >> ${out}/all_contigs.fa; done
if [[ ! -s ${out}/all_contigs.fa ]]; then error "something went wrong with joining files in $bin_folder into ${out}/all_contigs.fa"; fi


if [[ -s ${out}/megablast_out.raw.tab ]]; then 
	comm "megablast alignment already done. Skipping..."
else
	comm "aligning ${out}/all_contigs.fa to ${BLASTDB} database with MEGABLAST. This is the longest step - please be patient. You may look at the classification progress in ${out}/megablast_out.raw.tab"
	blastn -task megablast -num_threads $threads\
	 -db ${BLASTDB}/nt\
	 -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore evalue nident length'\
	 -query ${out}/all_contigs.fa > ${out}/megablast_out.raw.tab

	if [[ $? -ne 0 ]]; then error "Failed to run megablast. Exiting..."; fi
fi


comm "removing unnecessary lines that lead to bad tax IDs (without a proper rank)"
${SOFT}/prune_blast_hits.py ${TAXDUMP}/nodes.dmp ${out}/megablast_out.raw.tab > ${out}/megablast_out.pruned.tab
if [[ $? -ne 0 ]]; then error "Failed to run prune_blast_hits.py to remove the bad lines. Exiting..."; fi
cat ${out}/megablast_out.pruned.tab | cut -f1,2,3,4,5,7,8,9,10,11,12 > ${out}/megablast_out.tab
if [[ $? -ne 0 ]]; then error "Failed to remove extra column. Exiting..."; fi


comm "making mapping file"
cat ${out}/megablast_out.pruned.tab | cut -f5,6 > ${out}/mapping.tax
if [[ $? -ne 0 ]]; then error "Failed to make mapping file. Exiting..."; fi



########################################################################################################
########################             GET TAXONOMY WITH TAXATOR-TK               ########################
########################################################################################################
announcement "GET TAXONOMY FROM MEGABLAST OUTPUT WITH TAXATOR-TK"

export TAXATORTK_TAXONOMY_NCBI=$TAXDUMP

comm "pulling out classifications with taxator"
cat ${out}/megablast_out.tab | taxator -a megan-lca -t 0.3 -e 0.01 -g ${out}/mapping.tax > ${out}/predictions.gff3
if [[ $? -ne 0 ]]; then error "Failed to run taxator. Exiting..."; fi


comm "binning and consolidating classifications for each contig"
sort -k1,1 ${out}/predictions.gff3 | binner -n classification -i genus:0.6 > ${out}/binned_predictions.txt
if [[ $? -ne 0 ]]; then error "Failed to run binner. Exiting..."; fi
mv binning.log $out


comm "pulling out full taxonomy path with taxknife"
cat ${out}/binned_predictions.txt | taxknife -f 2 --mode annotate -s path | grep -v "Could not" | cut -f1,2 > ${out}/contig_taxonomy.tab
if [[ $? -ne 0 ]]; then error "Failed to extract full taxonomy path with taxknife. Exiting..."; fi


comm "finding consensus taxonomy for each bin"
${SOFT}/classify_bins.py ${out}/contig_taxonomy.tab $bin_folder > ${out}/bin_taxonomy.tab
if [[ $? -ne 0 ]]; then error "Failed to get consensus of bin taxonomy"; fi
cat ${out}/bin_taxonomy.tab


comm "you will find the consensus taxonomy of each bin in ${out}/bin_taxonomy.tab"

if false; then
comm "renaming bins to their best taxonomy"
mkdir ${out}/renamed_bins
for bin in $(ls $bin_folder); do 
	tax=$(cat ${out}/bin_taxonomy.tab | awk -v bin=$bin '$1==bin' | cut -f2)
	best_tax=$(echo $tax | rev | cut -d';' -f1 | rev | cut -d'_' -f1)
	no_spaces=${best_tax// /_}
	final_name=${no_spaces//./}
	echo "$bin was renamed to ${final_name}.fa"
	cp ${bin_folder}/$bin ${out}/renamed_bins/${final_name}.fa
done

comm "cleaning up..."
mkdir ${out}/work_files
for file in $(ls $out | grep -v "contig_taxonomy.tab" | grep -v "bin_taxonomy.tab"); do
	mv ${out}/$file ${out}/work_files/
done

comm "you will find the bins taxonomy estimates in ${out}/bin_taxonomy.tab"
fi


########################################################################################################
########################    CLASSIFICATION PIPELINE FINISHED SUCCESSFULLY!!!    ########################
########################################################################################################
announcement "BIN CLASSIFICATION PIPELINE FINISHED SUCCESSFULLY!!!"

