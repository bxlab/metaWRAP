#!/usr/bin/env perl

=head1 NAME

bowtie2_extract_reads_mapped_to_specific_contigs.pl

=head1 SYNOPSIS

bowtie2_extract_reads_mapped_to_specific_contigs.pl -s samfile -id contigids1.txt -id contigids2.txt [-u]

=head1 DESCRIPTION

- Takes a .sam file as input (you can provide a bam file using BASH subprocesses or pipes)
- Takes multiple files with contigids
- Assumes the bowtie2 .bam/sam file has the full read stored inside it (i.e uses soft clipping)
- Assumes reads are interleaved/paired in sam file
- Reads are returned in pairs. If a /1 read maps to one file and /2 to another file, then put both in /1's file
- If the same contigid if present in two files, both mapped read output files will contain those reads
  (i.e. reads will be duplicated. so it is the user's responsibility to ensure contig ids are not duplicated) 
- -u will also print unmapped.reads (reads that were either not mapped to ANY contig, or reads that were not mapped to
  any of the contigid lists

Needs:

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2012.04.28

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use IO::File;

my $samfile;
my @contigidfiles;
my $unmapped = "";
my $output_prefix = "";

GetOptions (
    "samfile=s"    => \$samfile,
    "idfiles=s{,}" => \@contigidfiles,
    "unmapped"     => \$unmapped,
    "out=s"        => \$output_prefix,
);

if (not $samfile or not @contigidfiles) {
    print STDERR "bowtie2_extract_reads_mapped_to_specific_contigs.pl -s samfile -id contigids1.txt -id contigids2.txt [-u]\n";
    exit;
}

#-------------------------------

#-------------------------------
# load each contigidfile (that we want to extract reads for) into memory
# create a filehandle to write to for each contigidfile called output_prefix.contigidfile.readnums

my %include_ids; # two level hash 
my %mapped_reads_fh;

for my $contigidfile (@contigidfiles) {
    $mapped_reads_fh{$contigidfile} = IO::File->new;
    open $mapped_reads_fh{$contigidfile}, "| gzip >$output_prefix$contigidfile.mapped.reads.fq.gz" or die $!;
    open ID, "<$contigidfile" or die $!;
    while (<ID>) {
        /^>?(\S+)/ and $include_ids{$contigidfile}{$1} = 1;
    }
}
$mapped_reads_fh{unmapped} = IO::File->new;
if ($unmapped) {
    open $mapped_reads_fh{unmapped}, "| gzip >${output_prefix}unmapped.reads.fq.gz" or die $!;
}

#-------------------------------
# Now, go through each line of the sam file
# and each time you see a contig id, check which of the contigidfiles it belongs in,
# and print the read to that file

open SAM, "<$samfile" or die $!;
my ($fr, $fc, $rr, $rc); # forward read, forward contig, etc...
while (<SAM>) {
    next if /^@/ or /^#/;
    # get forward read, forward contig, reverse read, reverse contig;
    my @F=split(/\t/);
    die "Does not seem to be a valid SAM file" unless $F[1] =~ /^\d+$/;
    if ($F[1] & 16) {
        $fr = "\@$F[0]\n" . &revcomp($F[9]) . "\n\+\n" . reverse($F[10]) . "\n";
    }
    else {
        $fr = "\@$F[0]\n$F[9]\n\+\n$F[10]\n";
    }
    $fc = $F[2];
    $_ = <SAM>;
    @F=split(/\t/);
    if ($F[1] & 16) {
        $rr = "\@$F[0]\n" . &revcomp($F[9]) . "\n\+\n" . reverse($F[10]) . "\n";
    }
    else {
        $rr = "\@$F[0]\n$F[9]\n\+\n$F[10]\n";
    }
    $rc = $F[2];
    if ($fc eq "*" and $rc eq "*" and $unmapped) {
        print {$mapped_reads_fh{unmapped}} "$fr$rr";
        next;
    }
    $fc = $rc if $fc eq "*";
    my $found = 0;
    for my $contigidfile (@contigidfiles) {
        if ($include_ids{$contigidfile}{$fc}) {
            print {$mapped_reads_fh{$contigidfile}} "$fr$rr";
            $found = 1;
        }
    }
   if (not $found and $unmapped) { 
       print {$mapped_reads_fh{unmapped}} "$fr$rr";
   }
}
close SAM;

#############################################################################

sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}

sub revcomp {
    $_ = shift @_;
    tr/atgcATGC/tacgTACG/;
    return scalar reverse $_;
}
