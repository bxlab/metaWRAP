#!/usr/bin/env bash
# this script creates the file with configuration paths to the metaWRAP modules.

echo "# MANDITORY! Paths to custon pipelines and scripts of metaWRAP" > bin/config-metawrap
echo "PIPES=$(pwd)/pipelines" >> bin/config-metawrap
echo "SOFT=$(pwd)/scripts" >> bin/config-metawrap
echo "" >> bin/config-metawrap

echo "# OPTIONAL databases (see 'Databases' section of metaWRAP README for details)" >> bin/config-metawrap
echo "# path to kraken standard database" >> bin/config-metawrap
echo "KRAKEN_DB=$(pwd)/kraken_standard_database" >> bin/config-metawrap
echo "" >> bin/config-metawrap

echo "# path to indexed human genome (see bmtagger website for guide). This insludes files hg38.bitmask and hg38.srprism.*"  >> bin/config-metawrap
echo "BMTAGGER_DB=$(pwd)/bmtagger_db" >> bin/config-metawrap
echo "" >> bin/config-metawrap

echo "# paths to BLAST databases" >> bin/config-metawrap
echo "BLASTDB=$(pwd)/NCBI_nt" >> bin/config-metawrap
echo "TAXDUMP=$(pwd)/NCBI_tax" >> bin/config-metawrap

echo "All paths configured! Please loot at the metaWRAP/bin/config-metawrap file to make sure everything is correct"
