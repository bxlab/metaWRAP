**New: The metaWRAP manuscript has been published in [Microbiome](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1).**

**Current metaWRAP version = 1.0.4**. To update, run `conda install metawrap-mg=1.0.4` (consider backing up your `bin/config-metawrap`file before you update)

# MetaWRAP - a flexible pipeline for genome-resolved metagenomic data analysis

 MetaWRAP aims to be an **easy-to-use metagenomic wrapper suite** that accomplishes the core tasks of metagenomic analysis from start to finish: read quality control, assembly, visualization, taxonomic profiling, extracting draft genomes (binning), and functional annotation. Additionally, metaWRAP takes bin extraction and analysis to the next level (see module overview below). While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple approach before you delve deeper into parameterization of your analysis. MetaWRAP can be applied to a variety of environments, including gut, water, and soil microbiomes (see metaWRAP paper for benchmarks). Each individual module of metaWRAP is a standalone program, which means you can use only the modules you are interested in for your data.
 
 ![General walkthrough of metaWRAP modules](https://i.imgur.com/6GqRsm3.png)

## Metagenomic bin recovery improvements

 In addition to being a tool wrapper, MetaWRAP offers a **powerful hybrid approach** for extracting high-quality draft genomes (bins) from metagenomic data by using a variety of software (metaBAT2, CONCOCT, and MaxBin2, for example, since they are already wrapped into the Binning module) and utilizing their individual strengths and minimizing their weaknesses. MetaWRAP's [bin refinement module](https://i.imgur.com/JL665Qo.png) outperforms not only individual binning approaches, but also other bin consolidation programs (Binning_refiner, DAS_Tool) in both synthetic and real datasets. I emphasize that because this module is a standalone component, I encourage you to use your favorite binning softwares for the 3 intitial predictions (they do not have to come from metaBAT2, CONCOCT and MaxBin2). These predictions can also come from different parameters of the same software.

![Bin_refinement performance comparison in different microbiome types](https://i.imgur.com/KSk3l2B.jpg)


 MetaWRAP also includes a novel [bin reassembly module](https://i.imgur.com/GUSMXl8.png), which allows to drastically improve the quality of a set of bins by extracting the reads belonging to each bin, and **reassembling the bins** with a more permissive, non-metagenomic assembler. In addition to improving the N50 of the bins, this modestly increases the completion of the bins, and drastically reduces contamination. I recommend you run the reassembly on the final bins set from the Bin_refinement module, but this can be any bin set.
 

## OVERVIEW OF METAWRAP MODULES:

#### Metagemonic data pre-processing modules:
	1) Read_QC: read trimming and human read removal
	2) Assembly: metagenomic assembly and QC with metaSPAdes or MegaHit
	3) Kraken: taxonomy profiling and visualization or reads or contigs
	
#### Bin processing modules:
	1) Binning: initial bin extraction with MaxBin2, metaBAT2, and/or CONCOCT
	2) Bin_refinement: consolidate of multiple binning predicitons into a superior bin set
	3) Reassemble_bins: reassemble bins to improve completion and N50, and reduce contamination
	4) Quant_bins: estimate bin abundance across samples
	5) Blobology: visualize the community and extracted bins with blobplots
	6) Classify_bins: conservative but accurate taxonomy prediction for bins
	7) Annotate_bins: functionally annotate genes in a set of bins
#### For more details, please consult the [metaWRAP module descriptions](https://github.com/bxlab/metaWRAP/blob/master/Module_descriptions.md) and the [publication preprint](https://www.biorxiv.org/content/early/2018/03/06/277442).


##  SYSTEM REQUIREMENTS
 The resource requirements for this pipeline will vary greatly based on the amount of data being processed, but due to large memory requirements of many software used (KRAKEN and metaSPAdes to name a few), I recommend at 8+ cores and 64GB+ RAM. MetaWRAP officially supports only Linux x64 systems, but may be manually installed on others.

## INSTALLATION

#### Basic installation:
 To start, download [miniconda2](https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh) and install it:
 ``` bash
 wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
 bash Miniconda2-latest-Linux-x86_64.sh
 ```
 
 Then add channels to your conda environment, and install metaWRAP (supports Linux64):
 ``` bash
 # ORDER IS IMPORTANT!!!
 conda config --add channels defaults
 conda config --add channels conda-forge
 conda config --add channels bioconda
 conda config --add channels ursky

 conda install -c ursky metawrap-mg
 ```
 
 #### Better installation
 The conda installation of metaWRAP will install over 140 software dependancies, which may cause some conflicts with your currenly installed packages. If you already use conda, it is strongly recommended to [set up a conda custom environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) and install metaWRAP only in there. That way your current conda environment and metaWRAP's environment do not not conflict.
``` bash
conda create -n metawrap-env python=2.7
source activate metawrap-env

# ORDER IS IMPORTANT!!!
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky

conda install -c ursky metawrap-mg
```

#### Manual installation:
 You may want to manually install metaWRAP if you want better control over your environment, if you are installing on a system other than Linux64, or you just really dislike conda. In any case, you will need to manually install the [relevant prerequisite programs](https://github.com/bxlab/metaWRAP/blob/master/installation/dependancies.md). When you are ready, download or clone this ripository, carefully configure the `metaWRAP/bin/config-metawrap` file, and add the `metaWRAP/bin/` directory to to the `$PATH`. Thats it! 
 
 
## DATABASES

 In addition to the Conda installation, you will need to configure the paths to some databases that you downloaded onto your system. Use your favorite text editor to configure these paths in /some/path/miniconda2/bin/config-metawrap and make sure everything looks correct. If you are unsure where this config file is, run:
 ``` bash
 which config-metawrap
 ```

This is very important if you want to use any functions requiring databases, but depending on what you plan to do, the databases are not mandatory for metaWRAP (see Database section below). [Follow this guide for download and configuration instructions](https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md).

|    Database     | Size  |  Used in module |
|:---------------:|:---------------:|:-----:| 
|Checkm_DB	 |1.4GB| binning, bin_refinement, reassemble_bins |
|KRAKEN standard database|161GB |  kraken |
| NCBI_nt |71GB |  blobology, classify_bins |
| NCBI_tax |283MB |  blobology, classify_bins |
|Indexed hg38  	|  20GB |  read_qc |


## DETAILED PIPELINE WALKTHROUGH

  ![Detailed pipeline walkthrough](https://i.imgur.com/HDUPeXC.png)


## USAGE

Please look at the [MetaWRAP usage tutorial](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md) for detailed run instructions and examples.

Once all the dependencies are in place, running metaWRAP is relatively simple. The main metaWRAP script wraps around all of its individual modules, which you can call independently.
```
metaWRAP -h
	Usage: metawrap [module] --help
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
metawrap assembly -h

Usage: metawrap assembly [options] -1 reads_1.fastq -2 reads_2.fastq -o output_dir
Options:

	-1 STR          forward fastq reads
	-2 STR          reverse fastq reads
	-o STR          output directory
	-m INT          memory in GB (default=10)
	-t INT          number of threads (defualt=1)

	--use-megahit		assemble with megahit (default)
	--use-metaspades	assemble with metaspades instead of megahit
```

### Citing metaWRAP
If you found metaWRAP usefull in your research, please cite the publication: [MetaWRAP - a flexible pipeline for genome-resolved metagenomic data analysis](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1).

### Acknowledgements
Author of pipeline: [Gherman Uritskiy](https://github.com/ursky).

Principal Investigators: [James Taylor](http://bio.jhu.edu/directory/james-taylor/) and [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 

All feedback is welcome! For errors and bugs, please open a new Issue thread on this github page, and I will try to get things patched as quickly as possible. Please include the version of metaWRAP you are using (run `metawrap -v`), For general questions, suggestions and other feedback, you can contact me at guritsk1@jhu.edu. 

