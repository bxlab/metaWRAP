```
MetaBAT Copyright (c) 2014, The Regents of the University of California, 
through Lawrence Berkeley National Laboratory (subject to receipt of any 
required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, 
please contact Berkeley Lab's technology transfer department at  TTD@lbl.gov 
referring to "MetaBAT (2014-075)."

NOTICE.  This software was developed under funding from the U.S. Department of 
Energy.  As such, the U.S. Government has been granted for itself and others 
acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in 
the Software to reproduce, prepare derivative works, and perform publicly and 
display publicly.  Beginning five (5) years after the date permission to assert 
copyright is obtained from the U.S. Department of Energy, and subject to any 
subsequent five (5) year renewals, the U.S. Government is granted for itself 
and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide 
license in the Software to reproduce, prepare derivative works, distribute 
copies to the public, perform publicly and display publicly, and to permit 
others to do so.

```
NOTICE:
------------
***Check out MetaBAT 2 Release Candidate on development branch below. (Tutorial will be available soon.)***

***Check out a new tutorial [Best Binning Practices](https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices). ***

***Check out the MetaBAT with ensemble binning. It often produces better bins. ***

***Check out the MetaBAT paper [here](https://peerj.com/articles/1165). ***

***Check out the [benchmark](https://bitbucket.org/berkeleylab/metabat/wiki/Benchmark_MetaHIT) for alternative software and [a large scale example](https://bitbucket.org/berkeleylab/metabat/wiki/Example_Large_Data).***

***Check out a [Nature paper](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature19366.html), where MetaBAT was used to build global ocean virome .***

***Check out a [Science paper](http://www.sciencemag.org/content/350/6259/434.abstract) using MetaBAT to identify novel methane producing archaeal genomes.***

***Check out a [Nature paper](http://www.nature.com/nature/journal/v521/n7551/abs/nature14238.html) using MetaBAT to identify novel genomes from permafrost samples.***

***MetaBAT is in [the list of top genomics papers](https://peerj.com/collections/26-genomics-papers-oct-2015/)!!***

INSTALLATION:
-------------

Requirements:

* boost >= 1.55.0
* python >= 2.7
* scons >= 2.1.0
* g++ >= 4.9
* zlib >= 1.2.4
* binutils >= 2.2.2

(samtools 1.2 is downloaded and installed automatically)

```
#!bash
#clean up old files
rm -f master.tar.gz
rm -f dev.tar.gz
rm -rf berkeleylab-metabat-*

#stable release version
wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
tar xzvf master.tar.gz
cd berkeleylab-metabat-*

#latest development version 
wget https://bitbucket.org/berkeleylab/metabat/get/dev.tar.gz
tar xzvf dev.tar.gz
cd berkeleylab-metabat-*

#run the installation script
scons install PREFIX=$HOME [BOOST_ROOT=$BOOST_ROOT]
```

See INSTALL.md for Operating System specific installation instructions

Binary distributions:
```
#!bash
# Linux 64-bit
wget https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-linux-x64_v0.32.4.tar.gz
tar xzvf metabat-static-binary-linux-x64_v0.32.4.tar.gz

# Mac OSX El Capitan 64-bit 
wget https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-osx-elcapitan-x64_v0.32.4.tar.gz
tar xzvf metabat-static-binary-osx-elcapitan-x64_v0.32.4.tar.gz

cd metabat
```

For technical supports, please send an email to ddkang@lbl.gov or rsegan@lbl.gov

DESCRIPTION:
------------
MetaBAT, An Efficient Tool for Accurately Reconstructing Single Genomes from Complex Microbial Communities 

Grouping large genomic fragments assembled from shotgun metagenomic sequences to deconvolute complex microbial communities, or metagenome binning, enables the study of individual organisms and their interactions. Here we developed an automated metagenome binning software, called MetaBAT, which integrates empirical probabilistic distances of genome abundance and tetranucleotide frequency. Tested on both synthetic and real metagenome datasets, MetaBAT outperforms alternative methods in both accuracy and computational efficiency. Applying MetaBAT to an assembly from 1,704 human gut samples formed 1,634 genome bins (>200kb) in 3 hours, where 621 genome bins are >50% complete with <5% contamination from other species. Further analysis shows that the quality of these genome bins approaches manually curated genomes. 

MetaBAT 2 USAGE: running on command line
--------------------------------

The easy way:
> runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]

The slightly less easy way:

a) Generate a depth file from assembly BAM files

>jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

b) Run metabat

>metabat2 -i assembly.fasta -a depth.txt -o bins_dir/bin

To check MetaBAT options:

> metabat2 -h

```
Allowed options:
  -h [ --help ]                     produce help message
  -i [ --inFile ] arg               Contigs in (gzipped) fasta file format [Mandatory]
  -o [ --outFile ] arg              Base file name and path for each bin. The default output is fasta format.
                                    Use -l option to output only contig names [Mandatory].
  -a [ --abdFile ] arg              A file having mean and variance of base coverage depth (tab delimited; 
                                    the first column should be contig names, and the first row will be 
                                    considered as the header and be skipped) [Optional].
  -m [ --minContig ] arg (=2500)    Minimum size of a contig for binning (should be >=1500).
  --maxP arg (=95)                  Percentage of 'good' contigs considered for binning decided by connection
                                    among contigs. The greater, the more sensitive.
  --minS arg (=60)                  Minimum score of a edge for binning (should be between 1 and 99). The 
                                    greater, the more specific.
  --maxEdges arg (=200)             Maximum number of edges per node. The greater, the more sensitive.
  --noAdd                           Turning off additional binning for lost or small contigs.
  --cvExt                           When a coverage file without variance (from third party tools) is used 
                                    instead of abdFile from jgi_summarize_bam_contig_depths.
  -x [ --minCV ] arg (=1)           Minimum mean coverage of a contig in each library for binning.
  --minCVSum arg (=1)               Minimum total effective mean coverage of a contig (sum of depth over 
                                    minCV) for binning.
  -s [ --minClsSize ] arg (=200000) Minimum size of a bin as the output.
  -t [ --numThreads ] arg (=0)      Number of threads to use (0: use all cores).
  -l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences.
  --saveTNF arg                     File to save (or load if exists) TNF matrix for each contig in input.
  --saveCls                         Save cluster memberships as a matrix format
  --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
  --noBinOut                        No bin output. Usually combined with --saveCls to check only contig 
                                    memberships
  -d [ --debug ]                    Debug output
  -v [ --verbose ]                  Verbose output
```
Choice of Options:

* In MetaBAT 2, parameter optimization will be unnecessary, though we allowed a few parameters so that advanced users might play with them. Most of given parameters are used to decide the amount of data used for binning.
* You can decrease -m (--minContig) when the qualities of both assembly and formed bins with default value are very good.
* You can decrease --maxP and --maxEdges when the qualities of both assembly and formed bins are very bad.
* You can increase --minS when the qualities of both assembly and formed bins are very bad.
* Set --noAdd when added small or leftover contigs cause too much contamination.

MetaBAT USAGE: running on command line
--------------------------------
* MetaBAT will be left as the legacy version since it might outperform MetaBAT 2 in certain conditions where there are many samples available and assembly quality is good. And with proper parameter tuning, it might outperform MetaBAT 2 in other conditions as well.  

Run metabat

>metabat -i assembly.fasta -a depth.txt -o bin

To check MetaBAT options:

> metabat -h

```
Allowed options:
  -h [ --help ]                     produce help message
  -i [ --inFile ] arg               Contigs in (gzipped) fasta file format [Mandatory]
  -o [ --outFile ] arg              Base file name for each bin. The default output is fasta format. Use -l
                                    option to output only contig names [Mandatory]
  -a [ --abdFile ] arg              A file having mean and variance of base coverage depth (tab delimited;
                                    the first column should be contig names, and the first row will be
                                    considered as the header and be skipped) [Optional]
  --cvExt                           When a coverage file without variance (from third party tools) is used
                                    instead of abdFile from jgi_summarize_bam_contig_depths
  -p [ --pairFile ] arg             A file having paired reads mapping information. Use it to increase
                                    sensitivity. (tab delimited; should have 3 columns of contig index
                                    (ordered by), its mate contig index, and supporting mean read coverage.
                                    The first row will be considered as the header and be skipped) [Optional]
  --p1 arg (=0)                     Probability cutoff for bin seeding. It mainly controls the number of
                                    potential bins and their specificity. The higher, the more (specific)
                                    bins would be. (Percentage; Should be between 0 and 100)
  --p2 arg (=0)                     Probability cutoff for secondary neighbors. It supports p1 and better be
                                    close to p1. (Percentage; Should be between 0 and 100)
  --minProb arg (=0)                Minimum probability for binning consideration. It controls sensitivity.
                                    Usually it should be >= 75. (Percentage; Should be between 0 and 100)
  --minBinned arg (=0)              Minimum proportion of already binned neighbors for one's membership
                                    inference. It contorls specificity. Usually it would be <= 50
                                    (Percentage; Should be between 0 and 100)
  --verysensitive                   For greater sensitivity, especially in a simple community. It is the
                                    shortcut for --p1 90 --p2 85 --pB 20 --minProb 75 --minBinned 20
                                    --minCorr 90
  --sensitive                       For better sensitivity [default]. It is the shortcut for --p1 90 --p2 90
                                    --pB 20 --minProb 80 --minBinned 40 --minCorr 92
  --specific                        For better specificity. Different from --sensitive when using correlation
                                    binning or ensemble binning. It is the shortcut for --p1 90 --p2 90 --pB
                                    30 --minProb 80 --minBinned 40 --minCorr 96
  --veryspecific                    For greater specificity. No correlation binning for short contig
                                    recruiting. It is the shortcut for --p1 90 --p2 90 --pB 40 --minProb 80
                                    --minBinned 40
  --superspecific                   For the best specificity. It is the shortcut for --p1 95 --p2 90 --pB 50
                                    --minProb 80 --minBinned 20
  --minCorr arg (=0)                Minimum pearson correlation coefficient for binning missed contigs to
                                    increase sensitivity (Helpful when there are many samples). Should be
                                    very high (>=90) to reduce contamination. (Percentage; Should be between
                                    0 and 100; 0 disables)
  --minSamples arg (=10)            Minimum number of sample sizes for considering correlation based
                                    recruiting
  -x [ --minCV ] arg (=1)           Minimum mean coverage of a contig to consider for abundance distance
                                    calculation in each library
  --minCVSum arg (=2)               Minimum total mean coverage of a contig (sum of all libraries) to
                                    consider for abundance distance calculation
  -s [ --minClsSize ] arg (=200000) Minimum size of a bin to be considered as the output
  -m [ --minContig ] arg (=2500)    Minimum size of a contig to be considered for binning (should be >=1500;
                                    ideally >=2500). If # of samples >= minSamples, small contigs (>=1000)
                                    will be given a chance to be recruited to existing bins by default.
  --minContigByCorr arg (=1000)     Minimum size of a contig to be considered for recruiting by pearson
                                    correlation coefficients (activated only if # of samples >= minSamples;
                                    disabled when minContigByCorr > minContig)
  -t [ --numThreads ] arg (=0)      Number of threads to use (0: use all cores)
  --minShared arg (=50)             Percentage cutoff for merging fuzzy contigs
  --fuzzy                           Binning with fuzziness which assigns multiple memberships of a contig to
                                    bins (activated only with --pairFile at the moment)
  -l [ --onlyLabel ]                Output only sequence labels as a list in a column without sequences
  -S [ --sumLowCV ]                 If set, then every sample that falls below the minCV will be used in an
                                    aggregate sample
  -V [ --maxVarRatio ] arg (=0)     Ignore any contigs where variance / mean exceeds this ratio (0 disables)
  --saveTNF arg                     File to save (or load if exists) TNF matrix for each contig in input
  --saveDistance arg                File to save (or load if exists) distance graph at lowest probability
                                    cutoff
  --saveCls                         Save cluster memberships as a matrix format
  --unbinned                        Generate [outFile].unbinned.fa file for unbinned contigs
  --noBinOut                        No bin output. Usually combined with --saveCls to check only contig
                                    memberships
  -B [ --B ] arg (=0)               Number of bootstrapping for ensemble binning (Recommended to be >=20)
  --pB arg (=50)                    Proportion of shared membership in bootstrapping. Major control for
                                    sensitivity/specificity. The higher, the specific. (Percentage; Should be
                                    between 0 and 100)
  --seed arg (=0)                   For reproducibility in ensemble binning, though it might produce slightly
                                    different results. (0: use random seed)
  --keep                            Keep the intermediate files for later usage
  -d [ --debug ]                    Debug output
  -v [ --verbose ]                  Verbose output
```

Choice of Options:

* '-i' input file should be either fasta or gzipped fasta file. (since v0.32.3)
* If '-a [--abdFile]' option is not given, TNF only binning will be executed.
* -p option is for utilizing paired info from short reads. It may improve sensitivity.
* '--p1' and '--p2' should be both high to maintain great specificity. Usually p1 >= p2 performs better.
* --minProb mainly controls the scope and sensitivity of binning. A smaller number improves sensitivity. It should be < p1, p2.
* --minBinned mainly controls the specificity. A greater number improves specificity. Usually <= 50. 
* Use --verysensitive on simple community for more inclusive binning.
* --minCorr would include contigs which are closely correlated in abundance but somewhat different in absolute abundance. More effective in availability of many samples. 
* Recruiting by correlation would be activated only if # of samples >= minSamples and be disabled (for better specificity) when minContigByCorr > minContig.
* --veryspecific and --superspecific will not recruit small contigs by abundance correlation.
* Coverage in a sample less than the number given with '-x [--minCV]' option will be ignored.
* Bin size less than the number given with '-s [--minClsSize]' option will not be reported.
* Contigs smaller than the length cutoff given with '-m [--minContig]' option will not be used in binning. The cutoff should be >= 1500; ideally >=2500.
* Smaller contigs (>1000) will be given a chance to be recruited to existing bins when # of samples >= minSamples by default setting.
* Use '-l [--onlyLabel]' when it is not necessary to record sequences. In this case, only the labels will be reported. In this case, use --noBinOut to prevent producing individual bins output.
* If '-S [--sumLowCV]' is set, the coverages smallers than minCV will be aggregated. It may improve performance in certain cases.
* If any number greater than 0 is given by '-V [--maxVarRatio]' option, contigs having spurious coverage pattern will be ignored. 
* '--saveDistance' option saves a lot of computations when multiple binning attempts are executed with different parameter settings.
* '--unbinned' option generates a file for unbinned contigs.
* '-B' option is for ensemble binning. Recommended to be 20 or more. Should be >= 10 for reasonable results. It tends to generate reduced number of better quality bins at the cost of some additional computation.
* '--pB' option controls for sensitivity and specificity tradeoff in ensemble binning. The smaller, the sensitive. Range is between 0 to 100. The default is 50.
* Produced bins would be stochastic if ensemble binning was used. --seed would minimize the stochasticity but still there would be slight difference.

OUTPUT:
-------
Each discovered bin will be saved as a fasta format

NOTES:
-------
The proper settings on the read aligner should be set to evenly distribute ambiguously mapping reads (the default option for bowtie2, bwa, and bbmap).

jgi_summarize_bam_contig_depths USAGE:
-----
> jgi_summarize_bam_contig_depths
```
Usage: jgi_summarize_bam_contig_depths <options> sortedBam1 [ sortedBam2 ...]
where options include:
	--outputDepth       arg  The file to put the contig by bam depth matrix (default: STDOUT)
	--percentIdentity   arg  The minimum end-to-end % identity of qualifying reads (default: 97)
	--pairedContigs     arg  The file to output the sparse matrix of contigs which paired reads span (default: none)
	--unmappedFastq     arg  The prefix to output unmapped reads from each bam file suffixed by 'bamfile.bam.fastq.gz'
	--noIntraDepthVariance   Do not include variance from mean depth along the contig
	--showDepth              Output a .depth file per bam for each contig base
	--includeEdgeBases       When calculating depth & variance, include the 1-readlength edges (off by default)

Options to control shredding contigs that are under represented by the reads
	--referenceFasta    arg  The reference file.  (It must be the same fasta that bams used)
	--shredLength       arg  The maximum length of the shreds
	--shredDepth        arg  The depth to generate overlapping shreds
	--minContigLength   arg  The mimimum length of contig to include for mapping and shredding
	--minContigDepth    arg  The minimum depth along contig at which to break the contig

```

Example with real data 
--------------------------------
Description of Data:

* 2 libraries of next-gen sequencing data of a mock community.
* The community is composed of 25 known genomes.
* Data is available to download: http://portal.nersc.gov/dna/RD/Metagenome_RD/MetaBAT/Software/Mockup/
* Refer to [benchmark](https://bitbucket.org/berkeleylab/metabat/wiki/Home) and [benchmark II](https://bitbucket.org/berkeleylab/metabat/wiki/Benchmark_MetaHIT) for more comprehensive examples.

After downloading the assembly file and two bam files, run the following command line (using MetaBAT >= v0.32.4):

>runMetaBat.sh -B 20 assembly.fa *.bam

MetaBAT forms about 28 bins. In this example, since we know the true membership of each contig, we can calculate the completeness and precision of a bin correctly.

```R
#The following is R commands (tested on Linux)
library(plyr)
library(foreach)
options(width=150)

refs <- read.table("membership.txt", sep="\t", header=T)

files <- list.files(".", pattern="*.[0-9]+.fa$", full.name=F)
bins <- foreach(f=files) %do% {
	system(sprintf("grep '>' %s | sed 's/>//'", f), intern=TRUE)
}

res <- foreach(b=bins, .combine=cbind) %do% {
	ddply(refs[match(b, refs$contig),], .(reference), function(x) sum(x$size), .drop=F)
}
rownames(res) <- res[,1]; res <- res[, seq(2,ncol(res),2)]
colnames(res) <- sapply(strsplit(files, "\\."), function(x) x[length(x)-1])

genome.size <- ddply(refs, .(reference), function(x) c(size=sum(x$size)))
genome.size <- genome.size[match(rownames(res), genome.size$reference),] 

res.prec <- apply(res, 2, function(x) max(x)/sum(x))
res.comp <- apply(res, 2, function(x) max(x)/genome.size$size[which.max(x)])
res.ref <- apply(res, 2, function(x) rownames(res)[which.max(x)])
res2 <- cbind.data.frame(Bin=names(res.prec), Size=colSums(res), 'Compl.'=res.comp, 'Prec.'=res.prec, 'Ref.'=res.ref)
res2
summary(res2[,2:4]) #Median bin size was 3.72Mb with 97% median completeness (or 80% by mean) and 100% median precision (or 96% by mean).
length(unique(res2$Ref.)) # 24 out of 25 genomes were binned.
```

The output will look like the following:

```
>res2
   Bin    Size     Compl.     Prec.                                           Ref.
1    1 7860794 0.97635160 0.4832147        Natronobacterium gregoryi SP2, DSM 3393
10  10 4115670 0.81602285 1.0000000              Spirochaeta smaragdinae DSM 11293
11  11 3921774 0.86124825 1.0000000           Natronococcus occultus SP4, DSM 3396
12  12 3874303 0.98891653 1.0000000                    Hirschia baltica ATCC 49814
13  13 3839986 0.99090122 1.0000000         Coraliomargarita akajimensis DSM 45221
14  14 3730710 0.98314815 1.0000000            Clostridium thermocellum ATCC 27405
15  15 3704412 0.87120817 1.0000000        Thermobacillus composti KWC4, DSM 18247
16  16 3684573 0.86379200 1.0000000          Frateuria aurantia Kondo 67, DSM 6220
17  17 3546737 0.98625152 1.0000000                  Meiothermus silvanus DSM 9946
18  18 3236335 0.99686126 1.0000000             Clostridium perfringens ATCC 13124
19  19 3173259 0.99905801 1.0000000                Segniliparus rotundus DSM 44985
2    2 5619150 0.99822442 1.0000000    Echinicola vietnamensis KMM 6221, DSM 17526
20  20 2924873 0.89827162 1.0000000          Corynebacterium glutamicum ATCC 13032
21  21 2147757 0.98982458 1.0000000    Fervidobacterium pennivorans Ven5, DSM 9078
22  22 2125683 0.98312892 1.0000000                         Olsenella uli DSM 7084
23  23 1805558 0.96656297 1.0000000                  Streptococcus pyogenes M1 GAS
24  24  703485 0.23799984 0.7072745                  Salmonella bongori NCTC 12419
25  25  576231 0.11425057 1.0000000              Spirochaeta smaragdinae DSM 11293
26  26  351671 0.06972657 1.0000000              Spirochaeta smaragdinae DSM 11293
27  27  292949 0.06889611 1.0000000        Thermobacillus composti KWC4, DSM 18247
28  28  274427 0.08428058 1.0000000          Corynebacterium glutamicum ATCC 13032
3    3 5272697 0.92686085 0.7153358      Escherichia coli str. K-12 substr. MG1655
4    4 5090248 0.97611998 1.0000000     Desulfotomaculum gibsoniae Groll, DSM 7213
5    5 4940344 0.99375427 1.0000000      Desulfosporosinus meridiei S10, DSM 13257
6    6 4829038 0.96953761 1.0000000   Desulfosporosinus acidophilus SJ4, DSM 22704
7    7 4766169 0.98037872 1.0000000           Terriglobus roseus KBS 63, DSM 18391
8    8 4753789 0.91043046 1.0000000 Salmonella enterica subsp. arizonae serovar 62
9    9 4735197 0.98375082 1.0000000                      Pseudomonas stutzeri RCH2

>summary(res2[,2:4])
      Size             Compl.           Prec.       
 Min.   : 274427   Min.   :0.0689   Min.   :0.4832  
 1st Qu.:2142238   1st Qu.:0.8632   1st Qu.:1.0000  
 Median :3717561   Median :0.9728   Median :1.0000  
 Mean   :3424922   Mean   :0.8031   Mean   :0.9609  
 3rd Qu.:4756884   3rd Qu.:0.9869   3rd Qu.:1.0000  
 Max.   :7860794   Max.   :0.9991   Max.   :1.0000  
```
MetaBAT can be run to the same data with different parameter settings. 
For instance, if one wants genome bins to be more sensitive (less specific), the following command can be used:

>runMetaBat.sh -B 20 --verysensitive assembly.fa *.bam