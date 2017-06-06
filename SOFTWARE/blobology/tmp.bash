#!/bin/bash




#HOW TO RUN:
#../../blobology.bash 12 assembly.fa reads_1.fastq reads_2.fastq




## Blobology pipeline. 2013-09-12 Sujai Kumar 
## sujai.kumar@gmail.com github.com/blaxterlab/blobology

## Needs:
## Read files (unless you are providing assembly as fasta and
## read alignments as BAM).
## ABySS in path (unless you are providing your own assembly)
## NCBI blast+ suite (2.2.28 or above, which provides taxids in hit)
## NCBI formatted nt blast database (because it includes taxon information)
## NCBI taxonomy dump (to calculate higher taxon levels from hit taxids)
## samtools in path (to manipulate bam files)
## github blaxterlab/blobology scripts in path
## fastq-mcf from the ea-utils suite (version 1.1.2-537) if you want to
## quality- and adapter- trim your raw reads

####==========================================================================
#### Download read files from SRA
####==========================================================================

## Uncomment the lines below if you are replicating the results in the
## Blobology paper. Otherwise provide your own reads/assembly in the assembly/
## alignment steps.
## ERR138445 is a 300 bp paired-end library
## ERR138446 is a 600 bp paired-end library

# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138445/ERR138445_1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138445/ERR138445_2.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138446/ERR138446_1.fastq.gz
# wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR138/ERR138446/ERR138446_2.fastq.gz

## Uncomment if you want to do quality- and adapter- trimming using fastq-mcf

#fastq-mcf -o >(gzip -c >ERR138445_1.mcf.fastq.gz) -o >(gzip -c >ERR138445_2.mcf.fastq.gz) -l 50 -q 20 --qual-mean 20 -R adapters.fa \
#    ERR138445_1.fastq.gz ERR138445_2.fastq.gz &>ERR138445.mcf.err
#fastq-mcf -o >(gzip -c >ERR138446_1.mcf.fastq.gz) -o >(gzip -c >ERR138446_2.mcf.fastq.gz) -l 50 -q 20 --qual-mean 20 -R adapters.fa \
#    ERR138446_1.fastq.gz ERR138446_2.fastq.gz &>ERR138446.mcf.err

####==========================================================================
#### GENERAL CONFIG VARIABLES
####==========================================================================

## Num of processors. Or, uncomment the line below and enter manually.

#NUMPROC=`grep "^processor" /proc/cpuinfo | tail -n 1 | awk '{print $3}'`
#NUMPROC=$1

## K-mer size for prelim assembly. Not needed if you are providing your own
## assembly fasta file.
## You can crudely estimate this for your data using
## http://dna.med.monash.edu.au/~torsten/velvet_advisor/

KMER=61

## Location of local installation of nt blast database
## (not needed if using blast remotely, which is slower).
## The NCBI nt databases can be downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/blast/db/ using the following command:

# wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
# for a in nt.*.tar.gz; do tar xzf $a; done

BLASTDB=/scratch/gu/RefSeq/NCBI_nt

## Location of NCBI tar gunzipped directory downloaded from
## ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
## Default: current directory

TAXDUMP=/scratch/gu/RefSeq/NCBI_tax

####==========================================================================
#### STEP 1, run ABySS assembler, or provide assembly fasta file
####==========================================================================

## If providing own assembly fasta file, uncomment and insert filename,
## and comment out the remaining commands in this section:

ASSEMBLY=$2

## If using ABySS, uncomment the commands below and insert your ownread files
## in the section 'lib="..."'.
## See http://www.bcgsc.ca/downloads/abyss/doc/#assemblingapaired-endlibrary
## for details. We also recommend running the
## command below on a cluster if possible (with np= the number of cores
## that you can request and j= number of processors on a single machine.
## e.g, if you have access to 6 nodes each with 8 processors, use np=48 j=8
## Even if you have only one machine with many cores, compile ABySS with
## mpi support so that it can run many processes in parallel.

#abyss-pe name=pa${KMER} k=${KMER} np=60 j=12 lib="lib1 lib2" lib1="ERR138445_1.mcf.fastq.gz ERR138445_2.mcf.fastq.gz" lib2="ERR138446_1.mcf.fastq.gz ERR138446_2.mcf.fastq.gz"

## At the end of the abyss step, the assembled sequence will be in
## ${NAME}-scaffolds.fa
## ABySS does not automatically discard short contigs, so keep only contigs >=200 bp:

#fastaqual_select.pl -l 200 -f ${NAME}-scaffolds.fa >${NAME}-scaffolds.fa.l200 

## Now, set the ASSEMBLY environment variable to that filename, so it can be used for subsequent steps

#ASSEMBLY=${NAME}-scaffolds.fa.l200

####==========================================================================
#### STEP 2, find best blast hits of a random sample of contigs
####==========================================================================

## the following command selects 10000 contigs at random:

RND_ASSEMBLY=random10k.$ASSEMBLY
#/home/guritsk1/programs/blobology/fastaqual_select.pl -f $ASSEMBLY -s r -n 10000 >$RND_ASSEMBLY
#cp $ASSEMBLY $RND_ASSEMBLY
## Change $BLASTDB to match your own path location of the nt blast database
## If running blast remotely on NCBI's servers, use -remote -db nt
## instead of -db $BLASTDB/nt
## -outfmt 6 qseqid staxids gives only a two-column output table with the
## query id and the taxid of the best hit

#blastn -task megablast -query random10k.$ASSEMBLY -db ${BLASTDB}/nt -evalue 1e-5 -num_threads $1 -max_target_seqs 1 -outfmt '6 qseqid sseqid staxids' -out tmp
#cut -f 1,3 tmp > random10k.assembly.nt.1e-5.megablast
#rm tmp


####==========================================================================
#### STEP 3, map reads back to assembly using bowtie2
####==========================================================================

## This step is not needed if you already have BAM files with the read alignments.
## The BAM files should have interleaved reads (i.e. read 1 followed by read 2, ...)
## even if the reads are mapped unpaired, because a later script will use that
## ordering to pick out read pairs even if only one read matched a given contig.

## By mapping the two libraries separately, we can check if the blobplots (TAGC plots)
## are library specific, which is a useful QC feature

#bowtie2-build $ASSEMBLY $ASSEMBLY

#/home/guritsk1/programs/blobology/shuffleSequences_fastx.pl 4 <(cat $3) <(cat $4) > tmp
#bowtie2 -x $ASSEMBLY --very-fast-local -k 1 -t -p $1 --reorder --mm -U tmp | samtools view -S -b -T $1 > ${ASSEMBLY%.*}.bowtie2.bam
#rm tmp

# Note: shuffleSequences_fastx.pl works on gunzipped files, hence the need to zcat inline

####==========================================================================
#### STEP 4, combine BAM files and blast hits to create data file for plotting
####==========================================================================

## ABySS will create one or more BAM files. Replace *.bam with own BAM files
## if you have aligned them separately without using ABySS

/home/guritsk1/programs/blobology/gc_cov_annotate.pl --blasttaxid random10k.assembly.nt.1e-5.megablast --assembly $ASSEMBLY --bam ${ASSEMBLY%.*}.bowtie2.bam --out blobplot.txt --taxdump $TAXDUMP --taxlist species order phylum superkingdom

# Note: change --taxdump . to point to the directory where you unpacked NCBI's taxdump.tar.gz
# Note: change --taxlist to any combination of NCBI's taxonomy levels (lowercase) such as class genus family, etc.
# The default is species order phylum superkingdom

####==========================================================================
#### STEP 4, create blobplots using R, or visualise using blobsplorer
####==========================================================================

## Tested with R 2.15.2 and ggplot2 0.9.3.1
## If you used a different --out option for naming the output data file in
## gc_cov_annotate.pl, use that in place of blobplot.txt below
## 0.01 is the threshold - if fewer than this proportion of annotated contigs
## have a particular annotation, they are not shown in the legend

#../../makeblobplot.R blobplot.txt 0.01 taxlevel_phylum

## The output file blobplot.txt can also be loaded into Blobsplorer - see github.com/mojones/blobsplorer

