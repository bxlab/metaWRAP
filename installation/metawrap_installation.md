#### To update to the latest version:
MetaWRAP is being constantly improved week to week as more bugs and issues pop up. Because of the scale of the project it is almost impossible to get a perfect working version as the dependancy software are constantly changing. I recommend to update to the newest version of metaWRAP on a monthly basis.

Before updating, back up your `config-metawrap` file so you do not have to re-do the database configurations. Then update with conda:
```
conda update -y -c ursky metawrap-mg

# or for a specific version:
conda install -y -c ursky metawrap-mg=1.2.2
```

If you are using the (recommended) manual instalation of metaWRAP, simply run `git pull` inside the metaWRAP directory.
It should also be noted that it is possible for th eupdates to produce strange behavior in complex conda environments, so if you experience issues the safest way is to just delete the old metawrap-env environment (`rm -r miniconda/envs/metawrap-env`) and re-install from scratch.  


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
conda install biopython blas=2.5 blast=2.6.0 bmtagger bowtie2 bwa checkm-genome fastqc kraken=1.1 krona=2.7 matplotlib maxbin2 megahit metabat2 pandas prokka quast r-ggplot2 r-recommended salmon samtools=1.9 seaborn spades trim-galore
# Note: this last solution is universal, but you may need to manually install concoct=1.0 and pplacer.
```

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
 # Note: may take a while
 
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
# Note: may take a while

 # To fix the CONCOCT endless warning messages in metaWRAP=1.2, run
 conda install -y blas=2.5=mkl
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
