###########################################################################
## Step 1. Create contaminant databases
###########################################################################

## After identifying the contaminants as being Bacterial, we can create a subset of the NCBI database that is limited to Bacteria only
## By doing a more specific blast against this database, and picking only those contigs that hit this database better than a nematode database,
## we can conservatively remove Bacterial contigs (and the reads that map to them)

# Create Bacteria database (taxon id 2):

curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=nuccore&dopt=gilist&qty=20000000&filter=all&term=txid2\[Organism:exp\]" >txid2.gids
blastdb_aliastool -dbtype nucl -gilist txid2.gids -db nt -out nt_Bacteria

# Create nematoda database (taxon id 6231):

curl "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&db=nuccore&dopt=gilist&qty=2000000&filter=all&term=txid6231\[Organism:exp\]" >txid6231.gids
blastdb_aliastool -dbtype nucl -gilist txid6231.gids -db nt -out nt_Nematoda

# Do the blasts:

blastn -task megablast -query $ASSEMBLY -db $BLASTDB/nt_Bacteria -outfmt 6 -evalue 1e-10 -max_target_seqs 1 -out $ASSEMBLY.nt_Bacteria.1e-10
blastn -task megablast -query $ASSEMBLY -db $BLASTDB/nt_Nematoda -outfmt 6 -evalue 1e-10 -max_target_seqs 1 -out $ASSEMBLY.nt_Nematoda.1e-10

# Separate by best blast hit:

blast_separate_taxa.pl -b1 $ASSEMBLY.nt_Nematoda.1e-10 -b2 $ASSEMBLY.nt_Bacteria.1e-10

# Make list of contigids to be removed:

# based on megablast 1e-10 hit:

cut -f1 $ASSEMBLY.nt_Bacteria.1e-10.only > toremove.contigids

# based on gc (col 3) and total cov (sum of cols 4 and 5):

awk '$3>=.55 && ($4+$5)<=50'  blobplot.txt | cut -f1 >>toremove.contigids
awk '$3>=.60 && ($4+$5)<=100' blobplot.txt | cut -f1 >>toremove.contigids
awk '$3>=.65'                 blobplot.txt | cut -f1 >>toremove.contigids
sort toremove.contigids | uniq | grep -v seqid >toremove.contigids.uniq; mv toremove.contigids.uniq toremove.contigids

# get list of contigs to keep by removing contaminant contigs from total list of all contig ids:
grep -v seqid blobplot.txt | cut -f1 | fgrep.pl - -v -f toremove.contigids >tokeep.contigids

# NOTE: fgrep.pl is like fgrep, but seems to work faster than the unix tool even though it's been written in perl

# Extract reads mapping to these contigs:

for LIBNAME in ERR138445 ERR138446
do
    bowtie2_extract_reads_mapped_to_specific_contigs.pl -s <(samtools view $ASSEMBLY.$LIBNAME.bowtie2.bam) -id tokeep.contigids -out $LIBNAME.
done
