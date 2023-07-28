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
	3) Kraken/Kraken2: taxonomy profiling and visualization or reads or contigs
	
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

#### Manual installation (this is best, if you are comfortable):
 The best way to install and manage metaWRAP is to install it directly from github, and then install all of its dependancies through conda. This is how I usually use metaWRAP, as it allows to easily update the versions of metawrap and other packages. This also works on MacOS as well as Unix.  
 
 0. Install mamba: `conda install -y mamba`. Mamba will effectively replace conda and do exactly the same thing, but _much_ faster.
 1. Download or clone this ripository: `git clone https://github.com/bxlab/metaWRAP.git`
 2. Carefully configure the `yourpath/metaWRAP/bin/config-metawrap` file to it points to your desired database locations (you can modify this later). Follow the [database configuration guide](https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md) for details.
 3. Make metaWRAP executable by adding `yourpath/metaWRAP/bin/` directory to to your `$PATH`. Either add the line `PATH=yourpath/metaWRAP/bin/:$PATH` to your `~/.bash_profile` script, or copy over the contents of `yourpath/metaWRAP/bin/` into a location already in your `$PATH` (such as `/usr/bin/` or `/miniconda2/bin/`). 
 4. (Optional but recommended) Make a new conda environment to install and manage all dependancies:
```
mamba create -y -n metawrap-env python=2.7
conda activate metawrap-env
```
5. Install all [metaWRAP dependancies](https://github.com/bxlab/metaWRAP/blob/master/conda_pkg/meta.yaml) with conda:
 ```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky

# Unix/Linux only
mamba install --only-deps -c ursky metawrap-mg
# `conda install --only-deps -c ursky metawrap-mg` also works, but much slower

# OR
mamba install biopython blas=2.5 blast=2.6.0 bmtagger bowtie2 bwa checkm-genome fastqc kraken=1.1 kraken=2.0 krona=2.7 matplotlib maxbin2 megahit metabat2 pandas prokka quast r-ggplot2 r-recommended salmon samtools=1.9 seaborn spades trim-galore minimap2
# Note: this last solution is more universal, but you may need to manually install concoct=1.0 and pplacer.
```

#### Express Conda/Mamba installation (the quickest but least configurable):
Directly create a metawrap-specific environment and install metawrap.
```
# install mamba (replaces conda, but much faster):
 conda install -y mamba 
 
# install metawrap:
 mamba create -y --name metawrap-env --channel ursky metawrap-mg=1.3.2
 conda activate metawrap-env

# To fix the CONCOCT endless warning messages in metaWRAP=1.2+, run
 conda install -y blas=2.5=mkl
```

#### Bioconda installation (avoid):
MetaWRAP is also available through the Bioconda channel. **However**, this distribution is not recommended for most users, as I will only push major releases to Bioconda (i.e. `v1.1`, `v1.2`). This source is meant for specific applications that require a Bioconda distribution. To get the latest version of metaWRAP with the newest patches and bug fixes, please install through the `-c ursky` channel, as seen above. 

```
# Bioconda installation (not recommended):
conda install -y mamba
mamba install -y -c bioconda metawrap

# or `conda install -y -c bioconda metawrap`
```

#### Docker installation (avoid):
If you are running on OSX and dont want to install manually, or prefer to work in containerized environments, then [Docker](https://quay.io/repository/biocontainers/metawrap?tab=info) could be the way to go. **However**, as with the Bioconda distribution, I will only push major releases to Bioconda (i.e. `v1.1`, `v1.2`). To get the latest version of metaWRAP with the newest patches and bug fixes, please install through the `-c ursky` channel, as seen above. If you still need to use Docker but run into bugs that have been fixed in the latest versions, you can manually update your scripts from this repository to apply the most recent patches. To install with Docker, run:
```
# Docker installation (not recommended and not fully supported)
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
| KRAKEN2 standard database|125GB | kraken2 |
| NCBI_nt |71GB |  blobology, classify_bins |
| NCBI_tax |283MB |  blobology, classify_bins |
| Indexed hg38  	|  20GB |  read_qc |


## DETAILED PIPELINE WALKTHROUGH

  ![Detailed pipeline walkthrough](https://i.imgur.com/HDUPeXC.png)
  Note: some features of this walkthrough are depricated since v0.7. To understand specific steps of each module, you can glance at the bash code in each script.


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
	kraken2		KRAKEN2 module
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
If you found metaWRAP usefull in your research, please cite the publication: [MetaWRAP - a flexible pipeline for genome-resolved metagenomic data analysis](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1). If certain software wrapped into metaWRAP were integral to your investigation (e.g. Salmon, MaxBin2, SPAdes, Kraken, etc.) please give them credit as well.

### Error reporting
The massive scale of the metaWRAP project unfortunately means that there are lots of opportunities for different components to fail depending on the exact environments it is installed on. Note that metaWRAP is simply a bash wrapper around other popular bioinformatics programs. If one of these other programs fails, the first thing to do is to troubleshoot the installation of that software, and not worry about metaWRAP itself until that component is fixed. If one of the components refuses to work on you environment, there may not be much I can do. Also remember that if you know a bit of bash/shell you can always see how metaWRAP calls these programs by investigating and possibly changing/tweaking the script files in `bin/metawrap-modules/`.

For errors and bugs relating to the actual metaWRAP software, please open a new Issue thread on this github page, however note that I no longer actively support metaWRAP due to moving on to other jobs. If you do report an error, please include the full output (stdout and stderr) from metaWRAP, and the version of you are using (run `metawrap -v`).

### Acknowledgements
Author of pipeline: [Gherman Uritskiy](https://github.com/ursky).

Principal Investigators: [James Taylor](http://bio.jhu.edu/directory/james-taylor/) and [Jocelyne DiRuggiero](http://bio.jhu.edu/directory/jocelyne-diruggiero/)

Institution: Johns Hopkins, [Department of Cell, Molecular, Developmental Biology, and Biophysics](http://cmdb.jhu.edu/) 
