# metaWRAP - Comprehensive Metagenome Analysis in One Place
## metaWRAP version 0.4

 metaWRAP aims to be an easy-to-use inclusive wrapper program that accomplishes the most basic tasks in metagenomic analysis: QC, assembly, binning, visualization, and taxonomic profiling. While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple first pass program before you delve deeper into parameterization of your approach. Each individual component of the pipeline is also a standalone module. This modularity allows the users to use only the modules they are interested in. 
 
 Additionally, metaWRAP offers a innovative hybrid pipeline for extracting high-quality draft genomes (bins) from metagenomic data. By using a variety of software (metaBAT2, CONCOCT, MaxBin2) and utilizing their individual strengths and minimizing their weaknesses, this pipeline will always produce stronger results than individual approaches. Additionally, due to its diverse binning approach, this pipeline shows promise to produce robust binning results in a variety of microbial communities. 

 If you already have your metagenomic data assembled and binned with two or more software, try using the BIN_REFINEMENT module to see how you can further improve your bin predictions! This modules uses Bin_refiner to make hybridized bin predictions based on two or three binning softwares, and then uses CheckM to consolidate all resulting bin sets to chose the best possible bins.
  

## OVERVIEW OF METAWRAP MODULES:

  There are 5 major parts to the pipeline:
  
    1) Read QC
    
    2) Assembly
    
    3) Binning super-module:
    	a. Binning with 3 softwares
		b. Bin refinement and consolidation
		c. Bin reassembly
		d. Bin quantitation
    
    4) Taxonomic profiling
    
    5) Visualization with Blobplots
  
  
    
  ![General walkthrough of metaWRAP modules](https://i.imgur.com/s9yAuQa.png)
   
    

  
## INSTALLATION

 Clone or download the metaWRAP directory into a semi-permanent location, then go into the metaWRAP/bin folder and edit the metaWRAP/bin/contig-metawrap file. Make sure that all the paths are correct, especially the paths pointing to the "scripts" and "pipelines" folders in the main metaWRAP directory. Once that is configured, simply copy the contents of metaWRAP/bin into your local bin folder, or simply add metaWRAP/bin/ to your path. If you are unsure how to do this, here are the commands to do this:
 
 ```
 echo "export PATH="/full/path/to/metaWRAP/bin:$PATH"" >> ~/.bash_profile_
 source ~/.bash_profile
 ```
 
 No try running metaWRAP -h to see if everything works!



## DEPENDENCIES

  Since this is a wrapper program, the biggest challenge in installing metaWRAP will likely be configuring all the dependencies correctly. Firstly, the path to the folder "meta-scripts", containing numerous scripts required for running this pipeline, needs to be configured in config.sh. Next, the following programs need to be installed in your PATH. NOTE: the versions of the programs may or may not be important. It is HIGHLY recommended that you to use miniconda to install everything. Look out for future releases or metaWRAP which will be directly wrapped into Conda. 
  
  NOTE: It is not necessary to install all of these, depending on which module you are interested in using. For example, if you do not want to sort out human reads, you dont have to install bmtagger and its database. Just make sure to use the --skip-bmtagger flag when running the read_qc module.

|    Software     | Tested version  |  Used in module 			|
|:---------------:|:---------------:|:---------------------------------:| 
|    BLAST        |    v=2.6.0      |  blobology			|
|    bmtagger     |    v=3.101      |  read_qc				|
|    Bowtie2      |    v=2.3.2      |  blobology			|
|    bwa          |    v=0.7.15     |  binning, reassemble_bins		|
|    Checkm       |    v=1.0.7      |  bin_refinement, reassemble_bins	|
|    FastQC       |    v=v0.11.5    |  read_qc				|
|    kraken       |    v=0.10.6     |  kraken				|
|    kronatools   |    v=2.7        |  kraken				|
|    megahit      |    v=1.1.1-2    |  assembly				|
|    SPAdes	  |    v=3.11.1	    |  assembly				|
|    metabat2     |    v=2.9.1      |  binning				|
|    concoct	  |    v=0.4.0	    |  binning				|
|    maxbin2      |    v=2.2.4      |  binning				|
|    perl         |    v=5.22.0     |  blobology			|
|    python       |    v=2.7.1      |  all modules			|
|    quast        |    v=4.5        |  assembly				|
|    R            |    v=3.3.2      |  blobology			|
|    samtools     |    v=1.3.1      |  assembly, blobology		|
|    SPAdes       |    v=3.10.1     |  assembly				|
|    trim_galore  |    v=0.4.3      |  read_qc     			|
|    salmon 	  |    v=0.8.1	    |  quant_bins			|


 To reiterate, the installation of most of these dependencies should not be difficult even in non-sudo environments with the use of conda. 


## DATABASES

Finally, you will need to download several databases and configure their paths in the config.sh file. This may be the longest step of the installation. Here is a full list of the databases:

|    Database     | Size  | Source |
|:---------------:|:---------------:|:-----:| 
|Checkm_DB		 |1.4GB| 	CheckM should prompt you to download this during first use	|
|KRAKEN standard database|161GB | 	look at the official KRAKEN support website for download instructions 		|
|RefSeq NCBI_nt 	|71GB | 	look at the config.sh for download instructions					|
|RefSeq NCBI_tax 	|283MB | 	look at the config.sh for download instructions					|
|Indexed hg38  		|  20GB | 	look at the bmtagger manual for instructions					|


  ![Detailed pipeline walkthrough](https://i.imgur.com/iNa6oUF.png)




## USAGE

Once all the dependencies are in place, running metaWRAP is relatively simple. The main metaWRAP script wraps around all of its individual modules, which you can call independently.

```
metaWRAP -h
	Usage: ./metaWRAP [module] --help
	Options:

	read_qc		Raw read QC module
	assembly	Assembly module
	binning		Binning module
	bin_refinement	Refinement of bins from binning module
	reassemble_bins Reassemble bins using metagenomic reads
	quant_bins	Quantify the abundance of each bin across samples
	blobology	Blobology module
	kraken		KRAKEN module
```

Each module is run separately. For example, to run the assembly module:

```
metaWRAP assembly -h

Usage: metaWRAP assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir
Options:

	-1 STR          forward fastq reads
	-2 STR          reverse fastq reads
	-o STR          output directory
	-m INT          memory in GB (default=10)
	-t INT          number of threads (defualt=1)

	--use-megahit		assemble with megahit (default)
	--use-metaspades	assemble with metaspades instead of megahit
```


Here is an example of running the binning module with multiple samples as input for better binning:
```bash
metaWRAP binning -t 48 -m 500 --checkm-best-bins --checkm-good-bins -a coassembly.fa -o binning_out sampleA_1.fastq sampleA_2.fastq sampleB_1.fastq sampleB_2.fastq sampleC_1.fastq sampleC_2.fastq
```



###  System requirements
 The resource requirements for this pipeline will vary greatly based on the number of reads you are processing, but I would advise against attempting to run it on anything less than 10 cores and 100GB RAM. With the help of conda, installing metaWRAP and its dependencies on a cluster (even without sudo privileges) should be relatively easy even for beginners.



### Acknowledgements
Author of pipeline: German Uritskiy.

Principal Investigators: [James Taylor](http://bio.jhu.edu/directory/james-taylor/) and [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

I do not claim to have any authorship of the many programs this pipeline uses. For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu, or leave a comment on this github page.

