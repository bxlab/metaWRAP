#!/usr/bin/env python2.7
# encoding:utf8
# authors: Erik Garrison, Sébastien Boisvert
"""This script takes two fastq or fastq.gz files and interleaves them
Usage:
    interleave-fasta fasta_file1 fasta_file2
"""

import sys

def interleave(f1, f2):
    """Interleaves two (open) fastq files.
    """
    while True:
        line = f1.readline()
        if line.strip() == "":
            break
        print line.strip()
        
        for i in xrange(3):
            print f1.readline().strip()
        
        for i in xrange(4):
            print f2.readline().strip()

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    except:
        print __doc__
        sys.exit(1)

    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1) as f1:
            with gzip.open(file2) as f2:
                interleave(f1, f2)
    else:
        with open(file1) as f1:
            with open(file2) as f2:
                interleave(f1, f2)

