# metaWRAP - Comprehensive Metagenomic Bin Analysis in One Place
## metaWRAP v=0.5

 MetaWRAP aims to be an easy-to-use inclusive wrapper program that accomplishes the most basic tasks in metagenomic analysis: QC, assembly, binning, visualization, and taxonomic profiling. While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple first pass program before you delve deeper into parameterization of your approach. Each individual component of the pipeline is also a standalone module. This modularity allows the users to use only the modules they are interested in. 
 
![General walkthrough of metaWRAP modules](https://i.imgur.com/LcC09ym.png)
   
 
 In addition to being a tool wrapper, MetaWRAP offers a innovative hybrid pipeline for extracting high-quality draft genomes (bins) from metagenomic data. By using a variety of software (metaBAT2, CONCOCT, MaxBin2) and utilizing their individual strengths and minimizing their weaknesses, the [bin refinement module](https://i.imgur.com/JL665Qo.png) will always produce stronger results than individual approaches. Additionally, due to its diverse binning approach, this pipeline shows promise to produce robust binning results in a variety of microbial communities. 

 MetaWRAP also includes a [bin reassembly module](https://i.imgur.com/GUSMXl8.png), which allows to drastically improve the quality of a set of bins by extracting the reads belonging to that draft genome, and reassembling it with a more permissive, non-metagenomic assembler. In addition to improving the N50 of the bins, this modestly increases the compleiton of the bins, and drastically reduces contamination.
  
 If you already have your metagenomic data assembled and binned with two or more software (or the same software with different parameters), try using the BIN_REFINEMENT and REASSEMBLE_BINS modules to see how you can further improve your bin predictions! 
  

## OVERVIEW OF METAWRAP MODULES:
  
#### Metagemonic data pre-processing modules:
		1) Read QC (trimming and human read removal)
    	2) Assembly (with metaSPAdes or MegaHit, plust assembly QC)
		3) Kraken (taxonomy profiling and visualization)
    	4) Binning (MaxBin2, metaBAT2, CONCOCT)
	
#### Bin processing modules:
		1) Bin refinement and consolidation of multiple bin sets
		2) Bin reassembly (reassemble bins to improve completiona and reduce contamination)
		3) Bin quantitation (bin abundance estimation across samples)
    	5) Blobology (visualize bin success with blobplots)
		6) Classify bins (asign taxonomy to draft genomes)
  

  
## INSTALLATION

 Clone or download the metaWRAP directory into a semi-permanent location, then go into the metaWRAP/bin folder and edit the metaWRAP/bin/contig-metawrap file. Make sure that all the paths are correct, especially the paths pointing to the "scripts" and "pipelines" folders in the main metaWRAP directory. Once that is configured, simply copy the contents of metaWRAP/bin into your local bin folder, or simply add metaWRAP/bin/ to your path:
 
 ```
 echo "export PATH="/full/path/to/metaWRAP/bin:$PATH"" >> ~/.bash_profile
 source ~/.bash_profile
 ```
 
 No try running ```metaWRAP -h``` to see if everything works!



## DEPENDENCIES

  Since this is a wrapper program, the biggest challenge in installing metaWRAP will likely be configuring all the dependencies correctly. A complete list of dependancies can be viewed in "dependancies.txt". Feel free to install them one by one. However, to make this as painless as possible, I HIGHLY recommend you use Conda, which automatically installs programs and handles dependancies. To start, download [miniconda2](https://conda.io/miniconda.html) (the Python 2.7 version) and install it. 
  
  Once you have miniconda2 installed, you can install most packages required by metaWRAP by using the metawrap-environment.yml file found in the metaWRAP folder. To create a new environment and install most of the dependancies, run:
  
  ``` bash
  conda env create -f metaWRAP/metawrap-environment.yml
  ```
  
  This will install over a hundred dependancy software for you in one go! They will be all added to the conda environment you just created. To use them (and to run metaWRAP), you need to enter the environment:
  
  ``` bash
  source activate metawrap
  ```
  When you are done with using metaWRAP and want to exit the environment, run:
  
  ``` bash
  source deactivate
  ```
  

  I want to emphasize that because metaWRAP is written in BASH scripts, it is realitively easy to understand why some part of the program is failing, or if there is a dependancy issue. Just to to that line in the code and see how metaWRAP is calling that software. There is not magic here - metaWRAP is simply calling programs from your current environment in a specific sequence just like you would if you were following a pipeline. If it cannot find them, an error pops up. Feel free to dive in the code to see where things went wrong. If you find a bug, please post in the "issues" section!


  After this, there are still a few software that you will need to install and place in yout PATH, becuase they are not currently in conda. Note that these software are only required only for specific modules. 

|    Software     | Tested version  |  Used in module 			|
|:---------------:|:---------------:|:---------------------------------:| 
|    kronatools   |    v=2.7        |  kraken				|
|    metabat2     |    v=2.9.1      |  binning				|
|    maxbin2      |    v=2.2.4      |  binning				|
|    pplacer      |    v=1.1        |  bin_refinement, reassemble_bins  |
|    taxator-kt   |    v=1.3.3      |  classify_bins                    |


## DATABASES

Finally, you will need to download several databases and configure their paths in the config.sh file. This may be the longest step of the installation. Again, you may not need all of these if you intend to use specific parts of the pipeline. Don't forget to configure the paths to them in the metaWRAP/bin/config-metawrap file! Here is a full list of the databases:

|    Database     | Size  |  Used in module |
|:---------------:|:---------------:|:-----:| 
|Checkm_DB	 |1.4GB| binning, bin_refinement, reassemble_bins |
|KRAKEN standard database|161GB |  kraken |
| NCBI_nt |71GB |  blobology |
| NCBI_tax |283MB |  blobology |
|Indexed hg38  	|  20GB |  read_qc |


### Downloading the CheckM database:
``` bash
mkdir MY_CHECKM_FOLDER
checkm data setRoot
# CheckM will prompt to to chose your storage location...
checkm data update

# If there is difficulty connecting to the servers, you can update manually:
cd MY_CHECKM_FOLDER
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_v1.0.9.tar.gz
tar -xvf checkm_data_v1.0.9.tar.gz
rm checkm_data_v1.0.9.tar.gz
```
Thats it! CheckM should know what folder to use as its database.

### Downloading the KRAKEN standard database:
Note: this will download the entire RefSeq database and index it, which takes a lot of computational power, storage space, and RAM. During database building, you will need >450GB of space and >250GB of RAM. With 24 cores, this will take >5 hours. Note that this is only needed if you intend on running the KRAKEN module.
``` bash
kraken-build --standard --threads 24 --db MY_KRAKEN_DATABASE
kraken-build --db MY_KRAKEN_DATABASE --clean
```
Do not forget to set the KRAKEN_DB variable in the contig-metawrap file in metaWRAP/bin/
``` bash
KRAKEN_DB=/path/to/my/database/MY_KRAKEN_DATABASE
```

### Downloading the NCBI_nt BLAST database:
``` bash
mkdir NCBI_nt
cd  NCBI_nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
for a in nt.*.tar.gz; do tar xzf $a; done
```
Do not forget to set the BLASTDB variable in the contig-metawrap file in metaWRAP/bin/
``` bash
BLASTDB=/your/location/of/database/NCBI_nt
```

### Downloading the NCBI_nt taxonomy:
``` bash
mkdir NCBI_tax
cd NCBI_tax
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```
Do not forget to set the TAXDUMP variable in the contig-metawrap file in metaWRAP/bin/
``` bash
TAXDUMP=/your/location/of/database/NCBI_tax
```

### Making human genome index for bmtagger
If you want to remove human reads from tour sequencing in the READ_QC module, you will need to dowlnoad and index the human genome. See the official bmtagger manual for detailed instructions: https://www.hmpdacc.org/hmp/doc/HumanSequenceRemoval_SOP.pdf

First, lets download and merge the human genome hg38:
``` bash 
mkdir BMTAGGER_INDEX
cd BMTAGGER_INDEX
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
gunzip *fa.gz
cat *fa > hg38.fa
rm chr*.fa
```
Now lets index the human genome. Note that the file names of the indeces must be exactly as specified for metaWRAP to recognize them! Also note that indexing will take considerable memory and time (here I pass 10GB of RAM as a -M parameter).
``` bash
bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 10000
```
Done! Now dont forget to specify the BMTAGGER_DB variable in the contig-metawrap file in metaWRAP/bin/
``` bash
BMTAGGER_DB=/path/to/your/index/BMTAGGER_INDEX
```

## DETAILED PIPELINE WALKTHROUGH

  ![Detailed pipeline walkthrough](https://i.imgur.com/5bb6vlY.jpg)



## USAGE

Once all the dependencies are in place, running metaWRAP is relatively simple. The main metaWRAP script wraps around all of its individual modules, which you can call independently.

```
metaWRAP -h
	Usage: metaWRAP [module] --help
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

