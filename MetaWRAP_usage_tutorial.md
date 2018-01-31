# Guide to analyzing metagenomic data with metaWRAP

## Step 0: Download sample metagenomic data from the metaHIT gut survey (or use your own unzipped, demultiplex, paired-end Illumina reads). 
```
# download data from 3 samples:
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_2.fastq.gz
```
```
# unzip the data
gunzip *qz
```
```
# place the raw sequencing reads into a new folder
mkdir RAW_READS
mv *fastq RAW_READS
ls RAW_READS
ERR011347_1.fastq
ERR011347_2.fastq
ERR011348_1.fastq
ERR011348_2.fastq
ERR011349_1.fastq
ERR011349_2.fastq
```


## Step 1: Run metaWRAP-Read_qc to trim the reads and remove human contamination
Note: you will need the bmtagger hg38 index to remove the human reads - see the metaWRAP database installation instructions. ALternatively, use the --skip-bmtagger flag of of the ReadQC module to only do the read trimming.
```
# individually process each sample
mkdir READ_QC
metawrap read_qc -1 RAW_READS/ERR011347_1.fastq -2 RAW_READS/ERR011347_2.fastq -t 24 -o READ_QC/ERR011347
metawrap read_qc -1 RAW_READS/ERR011348_1.fastq -2 RAW_READS/ERR011348_2.fastq -t 24 -o READ_QC/ERR011348
metawrap read_qc -1 RAW_READS/ERR011349_1.fastq -2 RAW_READS/ERR011349_2.fastq -t 24 -o READ_QC/ERR011349
```

Lets have a glance at one of the output folders: `ls READ_QC/ERR011347`
```
# these are html reports of the read quality before and after QC:
post-QC_report
pre-QC_report

# these are reads taht have been quality trimmed with trimgalore:
trimmed_1.fastq
trimmed_2.fastq

# these are human reads that have been pulled out with bmtagger (in case you need them):
human_reads_1.fastq
human_reads_2.fastq

# these are the final trimmed and de-contaminated reads:
final_pure_reads_1.fastq
final_pure_reads_2.fastq
```

```
# move over the final QC'ed reads into a new folder
mkdir CLEAN_READS
for i in READ_QC/*; do 
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done
```


## Step 2: Assembling the metagenomes with the metaWRAP-Assembly module
Note: Depending on your goals you may want to assemble each sample seperately, but for the purposes of analyzing the whole community across samples, we will be co-assembling our samples.
```
# concatinate the reads from all the samples:
cat CLEAN_READS/ERR*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/ERR*_2.fastq > CLEAN_READS/ALL_READS_2.fastq
```
```
# assemble the reads with metaSPAdes (usually prefered unless you have a very large data set):
metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 200 -t 96 --use-metaspades -o ASSEMBLY
```
You will find the assembly file in ASSEMBLY/final_assembly.fasta, and the QUAST assembly report html in ASSEMBLY/assembly_report.html!

```
# looking at the top 10 contigs shows we got some longer contigs (considering that we are working with just 7Gbp of data)!
grep ">" ASSEMBLY/final_assembly.fasta | head
>NODE_1_length_196124_cov_2.427049
>NODE_2_length_176373_cov_3.889994
>NODE_3_length_163601_cov_3.070200
>NODE_4_length_142996_cov_2.771017
>NODE_5_length_109931_cov_3.516837
>NODE_6_length_106321_cov_2.842875
>NODE_7_length_99368_cov_2.860703
>NODE_8_length_95669_cov_2.506714
>NODE_9_length_91511_cov_12.466716
>NODE_10_length_88949_cov_2.730882
```

## Step 3: Run Kraken module on both reads and the assembly
Running kraken on the reads will give us an idea of the taxonomic composition of the communities in the three samples, while running kraken on the assembly will give us an idea what taxonomic groups were assembled better than others (the assembly process is heavily biased and should not be used to infer overall community composition).
```
# run the Kraken module on all files at once, subsetting the reads to 1M reads per sample to speed up the run
metawrap kraken -o KRAKEN -t 96 -s 1000000 CLEAN_READS/ERR*fastq ASSEMBLY/final_assembly.fasta

# lets have a look out the output folder:
ERR011347.kraken  ERR011348.kraken  ERR011349.kraken  final_assembly.kraken  kronagram.html
ERR011347.krona   ERR011348.krona   ERR011349.krona   final_assembly.krona
```
The .kraken files contain the KRAKEN-estimated taxonomy of each read or contig, while the .krona files summarize taxonomy statistics to be fed into KronaTools, which makes the kronagram.html file. The kronagram.html contains all the taxonomy information from all the samples and the co-assembly. Inspecting the kronas in a web browser will show you what the community composition is like.



## Step 4: Bin the co-assembly with three different algorithms with the Binning module

The initial binning process with CONCOCT, MaxBin, and metaBAT will be the more time intensive steps (especially CONCOCT and MaxBin), so I would advise you to run the Binning module with each of the algorithms seperately. However, metaWRAP supports running all three together. Our dataset is reasonably small, so I will run all three binning predictions at the same time.

If you are used to a different binning software(s), feel free to run them instead. The downstream refinement process takes in up to 3 different bin sets.
```
# running the binning module with all three binners - notice how I put both the F and R read files at the end of the command.
metawrap binning -o INITIAL_BINNING -t 96 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/ERR*fastq
```













