# Guide to analyzing metagenomic data with metaWRAP

## Step 0: Download sample metagenomic data from the metaHIT gut survey (or use your own unzipped, demultiplex, paired-end Illumina reads). 

Download data from 3 samples:
```
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_2.fastq.gz
```

Unzip the data
```
gunzip *qz
```


Place the raw sequencing reads into a new folder
```
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
Note: you will need the bmtagger hg38 index to remove the human reads - see the metaWRAP database installation instructions. Alternatively, use the --skip-bmtagger flag of of the ReadQC module to only do the read trimming.

Individually process each sample
```
mkdir READ_QC
metawrap read_qc -1 RAW_READS/ERR011347_1.fastq -2 RAW_READS/ERR011347_2.fastq -t 24 -o READ_QC/ERR011347
metawrap read_qc -1 RAW_READS/ERR011348_1.fastq -2 RAW_READS/ERR011348_2.fastq -t 24 -o READ_QC/ERR011348
metawrap read_qc -1 RAW_READS/ERR011349_1.fastq -2 RAW_READS/ERR011349_2.fastq -t 24 -o READ_QC/ERR011349
```

Lets have a glance at one of the output folders: `ls READ_QC/ERR011347`

These are html reports of the read quality before and after QC:
```
post-QC_report
pre-QC_report
```

These are the final trimmed and de-contaminated reads:
```
final_pure_reads_1.fastq
final_pure_reads_2.fastq
```

Move over the final QC'ed reads into a new folder
```
mkdir CLEAN_READS
for i in READ_QC/*; do 
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done
```


## Step 2: Assembling the metagenomes with the metaWRAP-Assembly module
Note: Depending on your goals you may want to assemble each sample seperately, but for the purposes of analyzing the whole community across samples, we will be co-assembling our samples.

Concatinate the reads from all the samples:
```
cat CLEAN_READS/ERR*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/ERR*_2.fastq > CLEAN_READS/ALL_READS_2.fastq
```

Assemble the reads with the metaSPAdes option flag (usually prefered over MegaHIT unless you have a very large data set):
```
metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 200 -t 96 --use-metaspades -o ASSEMBLY
```

You will find the assembly file in ASSEMBLY/final_assembly.fasta, and the QUAST assembly report html in ASSEMBLY/assembly_report.html!

Looking at the top 10 contigs shows we got some longer contigs (considering that we are working with just 7Gbp of data)!
```
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

Run the Kraken module on all files at once, subsetting the reads to 1M reads per sample to speed up the run
```
metawrap kraken -o KRAKEN -t 96 -s 1000000 CLEAN_READS/ERR*fastq ASSEMBLY/final_assembly.fasta
```

Lets have a look out the output folder:
```
ERR011347.kraken  ERR011348.kraken  ERR011349.kraken  final_assembly.kraken  kronagram.html
ERR011347.krona   ERR011348.krona   ERR011349.krona   final_assembly.krona
```

The .kraken files contain the KRAKEN-estimated taxonomy of each read or contig, while the .krona files summarize taxonomy statistics to be fed into KronaTools, which makes the kronagram.html file. The kronagram.html contains all the taxonomy information from all the samples and the co-assembly. Inspecting the kronas in a web browser will show you what the community composition is like.



## Step 4: Bin the co-assembly with three different algorithms with the Binning module

The initial binning process with CONCOCT, MaxBin, and metaBAT will be the more time intensive steps (especially CONCOCT and MaxBin), so I would advise you to run the Binning module with each of the algorithms seperately. However, metaWRAP supports running all three together. Our dataset is reasonably small, so I will run all three binning predictions at the same time.

If you are used to a different binning software(s), feel free to run them instead. The downstream refinement process (the Bin_refinement module) takes in up to 3 different bin sets, although you can get around this by running splitting your bin sets into groups and then recursively consolidating them.

Run the binning module with all three binners - notice how I put both the F and R read files at the end of the command.
```
metawrap binning -o INITIAL_BINNING -t 96 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/ERR*fastq
```

In the output folder, we see folders with the 3 final bin sets, and a file containing the average and standard deviation of the library insert sizes (this may or may not be usefull to you).  
```
insert_sizes.txt  concoct_bins	maxbin2_bins  metabat2_bins  work_files
```
Looking inside these folders reveals that we found 47, 29, and 20 bins with concoct, metabat2, and maxbin2, respectively.


## Step 5: Consolidate 3 bin sets with the Bin_refinement module
Note: make sure you downloaded the CheckM database (see metaWRAP database instrucitons)

Now what you have metabat2, maxbin2, and concoct bins, lets consolidate them into a single, stronger bin set! If you used your own binning software, feel free to use any 3 bin sets. If you have more than 3, you can run them in groups. For example if you have 5 bin sets, try consolidating 1+2+3 and 4+5, and then consolidate again between the outputs.

When you do your refinement, make sure to put some thought into the minimum compleiton (-c) and maximum contamination (-x) parameters that you enter. During refinement, metaWRAP will have to chose the best version of each bin between 7 different versions of each bin. It will dynamically adjust to prioritize the bin quality that you desire. Consider this example: bin_123 comes in four versions in terms of completion/contamination: 95/15, 90/10, 80/5, 70/5. Which one is the best version? The high completion but high contamination, or the less complete but pure bin? This is subjective and depends on what you value in a bin and on what your purposes for bin extraction are. 

By default, the minimum completion if 70%, and maximum contamination is 5%. However, because of the relatively poor depth of these demonstration samples, we will set minimum completion to 50% and maximum contamination to 10%, but feel free to be much more picky. Parameters like -c 90 -x 5 are not unreasonable in some data (but you will get fewer bins, of course).

Run metaWRAP's Bin_refinement module:
```
metawrap bin_refinement -o BIN_REFINEMENT -t 96 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 50 -x 10
```







