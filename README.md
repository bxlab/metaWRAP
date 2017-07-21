## Introducing metaWRAP v0.1 - Comprehensive Metagenome Analysis for Beginners


 metaWRAP aims to be an easy-to-use inclusive wrapper program that accomplishes the most basic tasks in metagenomic analysis: QC, assembly, binning, visualization, and taxonomic profiling. While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple first pass program before you delve deeper into parameterization of your approach.
  
  There are 5 major parts to the pipeline:
  
    1) Read QC
    
    2) Assembly
    
    3) Binning
    
    4) Taxonomic profiling
    
    5) Visualization with Blobplots
    
    
  ![Detailed pipeline walkthrough](http://i.imgur.com/D1bOqLp.png)

  
## INSTALATION

 The core of metaWRAP is a bash script, so the program does not need any installation. Just download or clone the repository into your desired location. Just know that the core scripts need to be in the same location as the config.sh file, which contains all the paths to databases that metaWRAP uses. Once you have the scripts, go into the config.sh file, and edit all the paths to their correct locations. Have a look at the “DEPENDENCIES” section for more detail. Once everything is configured, you can run metaWRAP or any of the individual module scripts to display the help message.


## DEPENDENCIES

  Since this is a wrapper program, the biggest challenge in installing metaWRAP will likely be configuring all the dependencies correctly. Firstly, the path to the folder "meta-scripts", containing numerous scripts required for running this pipeline, needs to be configured in config.sh. Next, the following programs need to be installed in your PATH. NOTE: the versions of the programs may or may not be important.

|    Software     | Tested version  |
|:---------------:|:---------------:| 
|    BLAST        |    v=2.6.0      |
|    bmtagger     |    v=3.101      |
|    Bowtie2      |    v=2.3.2      |
|    bwa          |    v=0.7.15     |
|    Checkm       |    v=1.0.7      |
|    FastQC       |    v=v0.11.5    |
|    kraken       |    v=0.10.6     |
|    kronatools   |    v=2.7        |
|    megahit      |    v=1.1.1-2    |
|    metabat2     |    v=2.9.1      |
|    perl         |    v=5.22.0     |
|    python       |    v=2.7.1      |
|    quast        |    v=4.5        |
|    R            |    v=3.3.2      |
|    samtools     |    v=1.3.1      |
|    SPAdes       |    v=3.10.1     |
|    trim_galore  |    v=0.4.3      |

 The installation of most of these dependencies should not be difficult even in non-sudo environments with the use of conda. Future versions of metaWRAP are expected to have more detailed installation instructions, but for now you are on your own.


## DATABASES

Finally, you will need to download several databases and configure their paths in the config.sh file. This may be the longest step of the installation. Here is a full list of the databases:

|    Database     | Size  | Source |
|:---------------:|:---------------:|:-----:| 
|Checkm_DB		 |1.4GB| 	CheckM should prompt you to download this during first use	|
|KRAKEN standard database|161GB | 	look at the official KRAKEN support website for download instructions 		|
|RefSeq NCBI_nt 	|71GB | 	look at the config.sh for download instructions					|
|RefSeq NCBI_tax 	|283MB | 	look at the config.sh for download instructions					|
|Indexed hg38  		|  20GB | 	look at the bmtagger manual for instructions					|


## USAGE

Once all the dependencies are in place, running metaWRAP is relatively simple. You can chose to run metaWRAP itself, which will run all of the modules, or run individual modules as you wish. 

```
metaWRAP:
Usage: ./metaWRAP [options] -o output_folder raw_readsA_1.fastq raw_readsA_2.fastq
Options:
	-o STR          output directory
	-t INT          number of threads (default=1)
	-m INT          memory in GB (default=4)
	-1 STR		forward read file
	-2 STR		reverse read file
```


Here are examples of running metaWRAP and its sub-modules:
```bash
./metaWRAP –t 24 –m 120 -o metawrap_out -1 sample_1.fastq -2 sample_2.fastq 
./read_qc.sh –t 100 -1 reads_1.fastq -2 reads_2.fastq -o output_dir
./assembly.sh [options] -1 reads_1.fastq -2 reads_2.fastq -o assembly_dir -t 120 -m 800
./blobology.sh –t 10-m 20 -a final_assembly.fasta-1 reads_1.fastq -2 reads_2.fastq –o blobology_out
./binning.sh –t 48 –m 100 -a assembly.fa –o binning_out sampleA_1.fastq sampleA_2.fastq sampleX_1.fastq sampleX_2.fastq
./kraken.sh –t 24 –o kraken_out assembly.fasta reads_1.fastq reads_2.fastq
```


###  System requirements
 The resource requirements for this pipeline will vary greatly based on the number of reads you are processing, but I would advise against attempting to run it on anything less than 10 cores and 100GB RAM. With the help of conda, installing metaWRAP and its dependencies on a cluster (even without sudo privileges) should be relatively easy even for beginners.


### Acknowledgements
Author of pipeline: German Uritskiy.
Principal Investigators: [James Taylor](http://bio.jhu.edu/directory/james-taylor/) and [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)
Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

I do not claim to have any authorship of the many programs this pipeline uses. For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu, or leave a comment on this github page.


