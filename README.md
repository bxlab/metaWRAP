## Introducing metaWRAP v0.2 - Comprehensive Metagenome Analysis for Beginners


 metaWRAP aims to be an easy-to-use inclusive wrapper program that accomplishes the most basic tasks in metagenomic analysis: QC, assembly, binning, visualization, and taxonomic profiling. While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple first pass program before you delve deeper into parameterization of your approach.
  
  There are 5 major parts to the pipeline:
  
    1) Read QC
    
    2) Assembly
    
    3) Binning
    
    4) Taxonomic profiling
    
    5) Visualization with Blobplots
    
    
  ![Detailed pipeline walkthrough](http://i.imgur.com/D1bOqLp.png)

  
## INSTALATION

 Clone or download the metaWRAP directory into a semi-permanent location, then go into the metaWRAP/bin folder and edit the metaWRAP/bin/contig-metawrap file. Make sure that all the paths are correct, especially the paths pointing to the "scripts" and "pipelines" folders in the main metaWRAP directory. Once that is configured, simply copy the contents of the metaWRAP/bin into your local bin folder, or simply add them to your path. If youre unsure how to do the later, here are the commands:
 
 ```
 echo "export PATH="/full/path/to/metaWRAP/bin:$PATH"" >> ~/.bash_profile
 source ~/.bash_profile
 ```
 
 No try running metaWRAP -h to see if everything works!

## USAGE

Once all the dependencies are in place, running metaWRAP is relatively simple. The main metaWRAP script wraps around all of its indivirual modules, which you can call independantly.

```bash
metaWRAP -h

Usage: metaWRAP [module] --help
Options:

	all		Run comprehensive pipeline (all modules)
	read_qc		Raw read QC module
	assembly	Assembly module
	binning		Binning module
	blobology	Blobology module
	kraken		KRAKEN module
	phylosift	Phylosift module
```

To run all the modules automatically:

```bash
metaWRAP all -h

Usage: ./metaWRAP all [options] -o output_folder -1 raw_readsA_1.fastq -2 raw_readsA_2.fastq
Note: currently only works on one pair of reads. To run run multiple samples. However individual modules can handle more samples.
Options:

	-o STR          output directory
	-t INT          number of threads (default=1)
	-m INT          memory in GB (default=4)
	-1 STR		forward read file
	-2 STR		reverse read file
	--fast		skips some non-essential and computationally expensive parts of the pipeline:
				-skips human sequence removal and pre-qc fastq report in read_qc module
  	                        -replaces metaspades with megahit as primary assembler
          	                -subsamples Kraken and Phylosift modules to 10k reads instead of 1M reads
                  	        -skips bin reassembly in Binning module, and skips Kraken on bins
                          	-subsammples to 10k contigs and skips bin annotation in Blobology module
```


Or run each module seperately. For example, to run the assembly module:

```bash
metaWRAP assembly -h

Usage: metaWRAP assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir
Options:

	-1 STR          forward fastq reads
	-2 STR          reverse fastq reads
	-o STR          output directory
	-m INT          memory in GB (default=10)
	-t INT          number of threads (defualt=1)

	--metaspades-only	only assemble with metaspades
	--megahit-only		only assemble with megahit (much faster)
```


Here is an example of running metaWRAP and all of its sub-modules:
```bash
./metaWRAP all –t 24 –m 120 -o metawrap_out -1 sample_1.fastq -2 sample_2.fastq
```

Here is an example of running the binning module with multiple samples as input for better binning:
```bash
metaWRAP binning -t 48 -m 500 --checkm-best-bins --checkm-good-bins -a coassembly.fa -o binning_out sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq sampleC_1.fastq sampleC_2.fastq
```



## DEPENDENCIES

  Since this is a wrapper program, the biggest challenge in installing metaWRAP will likely be configuring all the dependencies correctly. Firstly, the path to the folder "meta-scripts", containing numerous scripts required for running this pipeline, needs to be configured in config.sh. Next, the following programs need to be installed in your PATH. NOTE: the versions of the programs may or may not be important.
  
  NOTE: It is not necessary to install all of these, depending on which module you are interested in using. For example, if you dont want to sort out human reads, you dont have to install bmtagger and its database. Just make sure to use the --skip-bmtagger flag when running the read_qc module.

|    Software     | Tested version  |  Used in module 		|
|:---------------:|:---------------:|:---------------------:| 
|    BLAST        |    v=2.6.0      |  blobology			|
|    bmtagger     |    v=3.101      |  read_qc				|
|    Bowtie2      |    v=2.3.2      |  blobology			|
|    bwa          |    v=0.7.15     |  binning				|
|    Checkm       |    v=1.0.7      |  binning				|
|    FastQC       |    v=v0.11.5    |  read_qc				|
|    kraken       |    v=0.10.6     |  kraken				|
|    kronatools   |    v=2.7        |  kraken				|
|    megahit      |    v=1.1.1-2    |  assembly				|
|    metabat2     |    v=2.9.1      |  binning				|
|    perl         |    v=5.22.0     |  blobology			|
|    python       |    v=2.7.1      |  all modules			|
|    quast        |    v=4.5        |  assembly				|
|    R            |    v=3.3.2      |  blobology			|
|    samtools     |    v=1.3.1      |  assembly, blobology	|
|    SPAdes       |    v=3.10.1     |  assembly				|
|    trim_galore  |    v=0.4.3      |  read_qc     			|

 The installation of most of these dependencies should not be difficult even in non-sudo environments with the use of conda. Future versions of metaWRAP are expected to have more detailed installation instructions, but for now you are on your own.


## DATABASES

Finally, you will need to download several databases and configure their paths in the config.sh file. This may be the longest step of the installation. However, you may skip databases that are not required for the modules/options you want to use. Here is a full list of the databases:

|    Database     | Size  | Source |
|:---------------:|:---------------:|:-----:| 
|Checkm_DB		 |1.4GB| 	CheckM should prompt you to download this during first use	|
|KRAKEN standard database|161GB | 	look at the official KRAKEN support website for download instructions 		|
|RefSeq NCBI_nt 	|71GB | 	look at the config.sh for download instructions					|
|RefSeq NCBI_tax 	|283MB | 	look at the config.sh for download instructions					|
|Indexed hg38  		|  20GB | 	look at the bmtagger manual for instructions					|



###  System requirements
 The resource requirements for this pipeline will vary greatly based on the number of reads you are processing, but I would advise against attempting to run it on anything less than 10 cores and 100GB RAM. With the help of conda, installing metaWRAP and its dependencies on a cluster (even without sudo privileges) should be relatively easy even for beginners.


### Acknowledgements
Author of pipeline: German Uritskiy.

Principal Investigators: [James Taylor](http://bio.jhu.edu/directory/james-taylor/) and [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

I do not claim to have any authorship of the many programs this pipeline uses. For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu, or leave a comment on this github page.

