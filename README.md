**To fix the CONCOCT endless warning messages in metaWRAP=1.2, run `conda install -y blas=2.5=mkl`

# MetaWRAP - a flexible pipeline for genome-resolved metagenomic data analysis

 MetaWRAP aims to be an **easy-to-use metagenomic wrapper suite** that accomplishes the core tasks of metagenomic analysis from start to finish: read quality control, assembly, visualization, taxonomic profiling, extracting draft genomes (binning), and functional annotation. Additionally, metaWRAP takes bin extraction and analysis to the next level (see module overview below). While there is no single best approach for processing metagenomic data, metaWRAP is meant to be a fast and simple approach before you delve deeper into parameterization of your analysis. MetaWRAP can be applied to a variety of environments, including gut, water, and soil microbiomes (see metaWRAP paper for benchmarks). Each individual module of metaWRAP is a standalone program, which means you can use only the modules you are interested in for your data.
 
 ![General walkthrough of metaWRAP modules](https://i.imgur.com/6GqRsm3.png)

## Metagenomic bin recovery improvements

 In addition to being a tool wrapper, MetaWRAP offers a **powerful hybrid approach** for extracting high-quality draft genomes (bins) from metagenomic data by using a variety of software (metaBAT2, CONCOCT, and MaxBin2, for example, since they are already wrapped into the Binning module) and utilizing their individual strengths and minimizing their weaknesses. MetaWRAP's [bin refinement module](https://i.imgur.com/JL665Qo.png) outperforms not only individual binning approaches, but also other bin consolidation programs (Binning_refiner, DAS_Tool) in both synthetic and real datasets. I emphasize that because this module is a standalone component, I encourage you to use your favorite binning softwares for the 3 intitial predictions (they do not have to come from metaBAT2, CONCOCT and MaxBin2). These predictions can also come from different parameters of the same software.

![Bin_refinement performance comparison in different microbiome types](https://i.imgur.com/KSk3l2B.jpg)


 MetaWRAP also includes a novel [bin reassembly module](https://i.imgur.com/GUSMXl8.png), which allows to drastically improve the quality of a set of bins by extracting the reads belonging to each bin, and **reassembling the bins** with a more permissive, non-metagenomic assembler. In addition to improving the N50 of the bins, this modestly increases the completion of the bins, and drastically reduces contamination. I recommend you run the reassembly on the final bins set from the Bin_refinement module, but this can be any bin set.
 

## OVERVIEW OF METAWRAP MODULES:

#### Metagemonic data pre-processing modules:
	1) Read_QC: read trimming and host (e.g. human) read removal
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
#### For more details, please consult the [metaWRAP module descriptions](https://github.com/bxlab/metaWRAP/blob/master/Module_descriptions.md) and the [publication](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1).


##  SYSTEM REQUIREMENTS
 The resource requirements for this pipeline will vary greatly based on the amount of data being processed, but due to large memory requirements of many software used (KRAKEN and metaSPAdes to name a few), I recommend at 8+ cores and 64GB+ RAM. MetaWRAP officially supports only Linux x64 systems, but may be installed on OSX manually or with docker (see below).

## INSTALLATION

#### To update to the latest version:
MetaWRAP is being constantly improved week to week as more bugs and issues pop up. Because of the scale of the project it is almost impossible to get a perfect working version as the dependancy software are constantly changing. I recommend to update to the newest version of metaWRAP on a monthly basis.

Before updating, back up your `config-metawrap` file so you do not have to re-do the database configurations. Then update with conda:
```
conda update -y -c ursky metawrap-mg

# or for a specific version:
conda install -y -c ursky metawrap-mg=1.2
```

If you are using the (recommended) manual instalation of metaWRAP, simply run `git pull` inside the metaWRAP directory.
It should also be noted that it is possible for th eupdates to produce strange behavior in complex conda environments, so if you experience issues the safest way is to just delete the old metawrap-env environment (`rm -r miniconda/envs/metawrap-env`) and re-install from scratch.  


#### Basic installation:
 To start, download [miniconda2](https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh) and install it:
 ``` bash
 wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
 bash Miniconda2-latest-Linux-x86_64.sh
 ```
 
 Then add channels to your conda environment, and install metaWRAP (supports Linux64):
 ``` bash
 # Note: ordering is important
 conda config --add channels defaults
 conda config --add channels conda-forge
 conda config --add channels bioconda
 conda config --add channels ursky

 conda install -y -c ursky metawrap-mg
 # Note: may take a couple hours
 
 # To fix the CONCOCT endless warning messages in metaWRAP=1.2, run
 conda install -y blas=2.5=mkl
 ```
 
 #### Better installation:
 The conda installation of metaWRAP will install over 140 software dependancies, which may cause some conflicts with your currenly installed packages. If you already use conda, it is strongly recommended to [set up a conda custom environment](https://conda.io/docs/user-guide/tasks/manage-environments.html) and install metaWRAP only in there. That way your current conda environment and metaWRAP's environment do not not conflict.
``` bash
conda create -y -n metawrap-env python=2.7
source activate metawrap-env

# Note: ordering is important
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky

conda install -y -c ursky metawrap-mg
# Note: may take a couple hours

 # To fix the CONCOCT endless warning messages in metaWRAP=1.2, run
 conda install -y blas=2.5=mkl
```

 
#### Best (manual) installation:
 This is how I usually use metaWRAP. By manually installing metaWRAP you will have better control over your environment and be able to change any programs and their versions as you see fit. You can also easily get the lastest updates throught `git pull`. If you are installing on a system other than Linux64, this may be your best bet. The hardest part is to install the [relevant prerequisite programs](https://github.com/bxlab/metaWRAP/blob/master/conda_pkg/meta.yaml), but this is much easier than you would think with the use of conda. Once you have these installed in your environment, download or clone this ripository, carefully configure the `metaWRAP/bin/config-metawrap` file, and add the `metaWRAP/bin/` directory to to the `$PATH` (likely by editing the `~/.bash_profile`). Alternatively, just copy over the `metaWRAP/bin/` contents into any location with excecutbale permission such as `/usr/bin/` or `/miniconda2/bin/` (depending on your permissions). Thats it!
 
 Easiest way to quickly install the dependancies:
 ```
 # Note: ordering is important
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky

# Unix/Linux only
conda install --only-deps -c ursky metawrap-mg

# OR
conda install biopython blas=2.5 blast=2.6.0 bmtagger bowtie2 bwa checkm-genome=1.0.12 concoct=1.0 fastqc kraken=1.1 krona=2.7 matplotlib maxbin2 megahit metabat2 pandas pplacer=1.1.alpha19 prokka quast r-ggplot2 r-recommended salmon samtools=1.9 seaborn spades trim-galore
 ```

#### Bioconda installation:
MetaWRAP is also available through the Bioconda channel. **However**, this distribution is not recommended for most users, as I will only push major releases to Bioconda (i.e. `v1.1`, `v1.2`). This source is meant for specific applications that require a Bioconda distribution. To get the latest version of metaWRAP with the newest patches and bug fixes, please install through the `-c ursky` channel, as seen above. 

```
# Bioconda installation (not recommended):
conda install -y -c bioconda metawrap
```

#### Docker installation:
If you are running on OSX and dont want to install manually, or prefer to work in containerized environments, then [Docker](https://quay.io/repository/biocontainers/metawrap?tab=info) could be the way to go. **However**, as with the Bioconda distribution, I will only push major releases to Bioconda (i.e. `v1.1`, `v1.2`). To get the latest version of metaWRAP with the newest patches and bug fixes, please install through the `-c ursky` channel, as seen above. If you still need to use Docker but run into bugs that have been fixed in the latest versions, you can manually update your scripts from this repository to apply the most recent patches. To install with Docker, run:
```
# Docker installation (not recommended)
docker pull quay.io/biocontainers/metawrap:1.2--0
```

 
## DATABASES

 In addition to the Conda installation, you will need to configure the paths to some databases that you downloaded onto your system. Use your favorite text editor to configure these paths in /some/path/miniconda2/bin/config-metawrap and make sure everything looks correct. If you are unsure where this config file is, run:
 ``` bash
 which config-metawrap
 ```

This is very important if you want to use any functions requiring databases, but depending on what you plan to do, the databases are not mandatory for metaWRAP (see Database section below). [Follow this guide for download and configuration instructions](https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md).

|    Database     | Size  |  Used in module |
|:---------------:|:---------------:|:-----:| 
| Checkm_DB |1.4GB| binning, bin_refinement, reassemble_bins |
| KRAKEN standard database|161GB |  kraken |
| NCBI_nt |71GB |  blobology, classify_bins |
| NCBI_tax |283MB |  blobology, classify_bins |
| Indexed hg38  	|  20GB |  read_qc |


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

All feedback is welcome! For errors and bugs, please open a new Issue thread on this github page, and I will try to get things patched as quickly as possible. Please include the **full** output (stdout and stderr) from metaWRAP, and the version of you are using (run `metawrap -v`), For general feedback you can contact me at guritsk1@jhu.edu, however please do not email me with error/bug reports. 

