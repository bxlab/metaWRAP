## ALL DEPENDANCIES:
- BLAST 2.6.0+
- bmtagger 3.101
- Bowtie2 2.3.2
- bwa 0.7.15-r1140
- Checkm v1.0.7
- FastQC v0.11.5
- kraken 0.10.6
- kronatools 2.7
- megahit v1.1.1-2-g02102e1
- metabat2 2.9.1
- concoct 0.4.0
- MaxBin2 2.2.4
- perl v5.22.0
- python2.7
- quast v4.5
- R 3.3.2
- ggplot2
- samtools 1.3.1
- SPAdes v3.10.1
- trim_galore 0.4.3
- python 2.7
- seaborn 0.8.1
- salmon
- taxator-tk
- prokka
- minimap2 2.24

## More detailed dependancies:

### general:
- python2.7
- shell

### assembly
- SPAdes v3.10.1
- bwa 0.7.15-r1140
- megahit v1.1.1-2-g02102e1
- quast v4.5

### binning
- bwa 0.7.15-r1140
- samtools 1.3.1
- metabat2 2.9.1
- concoct 0.4.0
- MaxBin2 2.2.4

#### bin_refinement
- bwa 0.7.15-r1140
- samtools 1.3.1
- kraken 0.10.6
- kraken_db STANDARD_DATABASE
- kronatools 2.7
- checkm v1.0.7
- checkm_DB (standard)
- SPAdes v3.10.1

#### bin_reassembly
- minimap2 2.24

### quant_bins
- salmon
- seaborn 0.8.1

### blobology
- refseq_db NCBI_nt
- refseq_db NCBI_tax
- blast 2.6.0+
- bowtie2 2.3.2
- perl v5.22.0
- R 3.3.2
- ggplot2
 
### kraken
- kraken 0.10.6
- kraken_db STANDARD_DATABASE
- kronatools 2.7

### read_qc
- fastqc v0.11.5
- trim_galore 0.4.3
- bmtagger 3.101
- human_db hg38

### classify_bins
- taxator-kt
- blastn
- NCBI_nt database
- NCBI_tax database

### annotate_bins
- prokka
