#!/bin/bash
# Usage: deinterleave_fastq.sh < interleaved.fastq f.fastq r.fastq [compress]
# 
# Deinterleaves a FASTQ file of paired reads into two FASTQ
# files specified on the command line. Optionally GZip compresses the output
# FASTQ files using pigz if the 3rd command line argument is the word "compress"
# 
# Can deinterleave 100 million paired reads (200 million total
# reads; a 43Gbyte file), in memory (/dev/shm), in 4m15s (255s)
# 
# Latest code: https://gist.github.com/3521724
# Also see my interleaving script: https://gist.github.com/4544979
# 
# Inspired by Torsten Seemann's blog post:
# http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

# Set up some defaults
GZIP_OUTPUT=0
PIGZ_COMPRESSION_THREADS=10

# If the third argument is the word "compress" then we'll compress the output using pigz
if [[ $3 == "compress" ]]; then
  GZIP_OUTPUT=1
fi

if [[ ${GZIP_OUTPUT} == 0 ]]; then
  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > $1) | cut -f 5-8 | tr "\t" "\n" > $2
else
  paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $1) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $2
fi
