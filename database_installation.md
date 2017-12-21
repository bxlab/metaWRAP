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
