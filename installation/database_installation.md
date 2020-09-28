## Downloading the CheckM database:
``` bash
mkdir MY_CHECKM_FOLDER

# Now manually download the database:
cd MY_CHECKM_FOLDER
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
cd ../

# Now you need to tell CheckM where to find this data before running anything:
checkm data setRoot     # CheckM will prompt to to chose your storage location
# On newer versions of CheckM, you would run:
checkm data setRoot /path/to/your/dir/MY_CHECKM_FOLDER
```
Thats it! CheckM should now use that folder and its contents as its database.



## Downloading the KRAKEN1 standard database:
Note: As of metaWRAP v1.3.2, we recomend you use Kraken2 instead of the original Kraken1 (see below).
Note: this will download the entire RefSeq database and index it, which takes a lot of computational power, storage space, and RAM. During database building, you will need >450GB of space and >250GB of RAM. With 24 cores, this will take >5 hours. Note that this is only needed if you intend on running the KRAKEN module.
``` bash
kraken-build --standard --threads 24 --db MY_KRAKEN_DATABASE
kraken-build --db MY_KRAKEN_DATABASE --clean
```
Do not forget to set the KRAKEN_DB variable in the config-metawrap file! Run `which config-metawrap` to find it.
``` bash
KRAKEN_DB=/path/to/my/database/MY_KRAKEN_DATABASE
```

## Downloading the KRAKEN2 standard database:
Note: Compared to Kraken1, the Kraken2 database is considerably more compact, making the download and indexing process much faster and less taxing on the system. You will need an estimated 120GB of RAM and ~128GB of space. Note that this is only needed if you intend on running the KRAKEN2 module.
``` bash
kraken2-build --standard --threads 24 --db MY_KRAKEN2_DB
```
Do not forget to set the KRAKEN_DB variable in the config-metawrap file! Run `which config-metawrap` to find it.
``` bash
KRAKEN2_DB=/path/to/my/database/MY_KRAKEN2_DATABASE
```

## Downloading the NCBI_nt BLAST database:
``` bash
mkdir NCBI_nt
cd  NCBI_nt
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz"
for a in nt.*.tar.gz; do tar xzf $a; done
```
Note: if you are using a more recent blast verions (beyond v2.6) you will need a the newer database format: `wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v4/nt_v4.*.tar.gz"`

Do not forget to set the BLASTDB variable in the config-metawrap file!
``` bash
BLASTDB=/your/location/of/database/NCBI_nt
```


## Downloading the NCBI taxonomy:
``` bash
mkdir NCBI_tax
cd NCBI_tax
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvf taxdump.tar.gz
```
Do not forget to set the TAXDUMP variable in the config-metawrap file! Run `which config-metawrap` to find it.
``` bash
TAXDUMP=/your/location/of/database/NCBI_tax
```


## Making host genome index for bmtagger
If you want to remove human (see end of instrutions for non-human hosts) reads from tour sequencing in the READ_QC module, you will need to dowlnoad and index the human genome. See the official bmtagger manual for detailed instructions: https://www.hmpdacc.org/hmp/doc/HumanSequenceRemoval_SOP.pdf

First, lets download and merge the human genome hg38:
``` bash 
mkdir BMTAGGER_INDEX
cd BMTAGGER_INDEX
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
gunzip *fa.gz
cat *fa > hg38.fa
rm chr*.fa
```
Now lets index the human genome. Note that the file names of the indeces must be exactly as specified for metaWRAP to recognize them! Also note that indexing will take considerable memory and time (here I pass 100GB of RAM as a -M parameter).
``` bash
bmtool -d hg38.fa -o hg38.bitmask
srprism mkindex -i hg38.fa -o hg38.srprism -M 100000
```
Note: metaWRAP looks for files hg38.bitmask and hg38.srprism - make sure they are names exactly like this.
Done! Now dont forget to specify the BMTAGGER_DB variable in the config-metawrap file! Run `which config-metawrap` to find it.
``` bash
BMTAGGER_DB=/path/to/your/index/BMTAGGER_INDEX
```

### Instructions for non-human hosts: 
For non-human hosts, the protocol for building and using the bmtagger index should be the same. In the end, you need to have th `.srprism` and `.bitmask` index files in the `BMTAGGER_INDEX` directory that you will link to in the `config-metawrap` file. When you run `metawrap read_qc`, you need to use the `-x` option to specify the prefix of your host. For example, if your `BMTAGGER_INDEX` directory has files `pig.srprism` and `pig.bitmask`, use `-x pig` as the option. 

