# MetaWRAP - Wrapper for Metagenomic Bin Analysis
## MetaWRAP v=0.6

 MetaWRAP aims to be an easy-to-use wrapper program that accomplishes the most basic tasks in metagenomic analysis: QC, assembly, binning, visualization, and taxonomic profiling. While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple first pass program before you delve deeper into parameterization of your approach. Each individual module of the pipeline is also a standalone component.
 
 In addition to being a tool wrapper, MetaWRAP offers a powerfull hybrid approach for extracting high-quality draft genomes (bins) from metagenomic data by using a variety of software (metaBAT2, CONCOCT, MaxBin2) and utilizing their individual strengths and minimizing their weaknesses. MetaWRAP's [bin refinement module](https://i.imgur.com/JL665Qo.png) outperforms not only individual binning approaches, but also other bin consolidation programs (Binning_refiner, DAS_Tool) in both synthetic and real datasets.

 MetaWRAP also includes a novel [bin reassembly module](https://i.imgur.com/GUSMXl8.png), which allows to drastically improve the quality of a set of bins by extracting the reads belonging to that draft genome, and reassembling it with a more permissive, non-metagenomic assembler. In addition to improving the N50 of the bins, this modestly increases the compleiton of the bins, and drastically reduces contamination.
 

## OVERVIEW OF METAWRAP MODULES:
![General walkthrough of metaWRAP modules](https://i.imgur.com/LcC09ym.png)

#### Metagemonic data pre-processing modules:
	1) Read QC (trimming and human read removal)
	2) Assembly (with metaSPAdes or MegaHit, plust assembly QC)
	3) Kraken (taxonomy profiling and visualization)
	4) Binning (MaxBin2, metaBAT2, CONCOCT)
	
#### Bin processing modules:
	1) Bin refinement and consolidation of multiple bin sets
	2) Bin reassembly (reassemble bins to improve completion and reduce contamination)
	3) Bin quantitation (bin abundance estimation across samples)
	5) Blobology (visualize bin success with blobplots)
	6) Classify bins (assign taxonomy to draft genomes)

##  SYSTEM REQUIREMENTS
 The resource requirements for this pipeline will vary greatly based on the amount of data being processed, but due to large memory requirements of many software used (KRAKEN and metaSPAdes to name a few), I would advise against attempting to run it on anything less than 10 cores and 100GB RAM. MetaWRAP officially supports only Linux x64 systems.


## INSTALLATION
 To start, download [miniconda2](https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh) and install it. Then add channels to your conda environment:
 ``` bash
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

 Once you have conda installed and configured, you can install metawrap and all its dependancies with the following command (unly supports Linux x64):
 ``` bash
 conda install -c ursky metawrap-binning
 ```
 WARNING: This will install over 140 software dependancies, which may cause some conflicts with your currenly installed packages. If you already use conda, it is recommended to [set up a conda custom environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) and install metaWRAP only in there. That way your current conda environment and metaWRAP's environment do not not conflict.
 ``` bash
 conda create -n metawrap-env python=2.7
 source activate metawrap-env
 
 conda config --add channels r
 conda config --add channels defaults
 conda config --add channels conda-forge
 conda config --add channels bioconda

 conda install -c ursky metawrap-binning
 ```
 
 If everything went well, running the following command should result in a help message
 ``` bash
 metaWRAP read_qc -h
 ```
 
 
## DATABASES

 Finally, use your favorite text editor to configure paths to databases in miniconda2/bin/config-metawrap and make sure all the paths look correct. This is very important if you want to use databases (see Database section below). If you are unsure where this config file is, run:
 ``` bash
 which config-metawrap
 ```

You will need to download and configure several databases and adjust their paths in the config-metawrap file. Note that depending on what modules you plan on using, you may not need all the databases. [Follow this guide for downloading and configuration instructions](https://github.com/ursky/metaWRAP/blob/master/installation/database_installation.md).

|    Database     | Size  |  Used in module |
|:---------------:|:---------------:|:-----:| 
|Checkm_DB	 |1.4GB| binning, bin_refinement, reassemble_bins |
|KRAKEN standard database|161GB |  kraken |
| NCBI_nt |71GB |  blobology, classify_bins |
| NCBI_tax |283MB |  blobology, classify_bins |
|Indexed hg38  	|  20GB |  read_qc |


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

### Acknowledgements
Author of pipeline: [German Uritskiy](https://github.com/ursky).

Principal Investigators: [James Taylor](http://bio.jhu.edu/directory/james-taylor/) and [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

All feedback is welcome! For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu, or leave a comment on this github page.

