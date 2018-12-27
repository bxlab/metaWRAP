
#### Basic installation:
 To start, download [miniconda2](https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh) and install it. Then add channels to your conda environment, and install metaWRAP (supports Linux64):
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
 
