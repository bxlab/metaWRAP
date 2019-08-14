#!/usr/bin/env perl

=head1 NAME
gc_cov_annotate.pl
=head1 SYNOPSIS
gc_cov_annotate.pl --blasttaxid CONTIGTAXIDFILE --assembly ASSEMBLYFASTAFILE [--taxdump TAXDUMPDIR] [--bam BAMFILE...]
=head1 AUTHORS
sujai.kumar@zoo.ox.ac.uk 2013.09.15
=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($blasttaxid_file, $taxdump_dir, $assembly_file, $output_file) = ("",".","","");
my @tax_list;
my @bam_files;
my @cov_files;

GetOptions (
    "blasttaxid=s" => \$blasttaxid_file,
    "assembly=s"   => \$assembly_file,
    "out:s"        => \$output_file,
    "taxdump:s"    => \$taxdump_dir,
    "bam:s{,}"     => \@bam_files,
    "cov:s{,}"     => \@cov_files,
    "taxlist:s{,}" => \@tax_list,
);

if (not @tax_list) {@tax_list = ("species","order","phylum","superkingdom")};
my %tax_levels;
foreach (@tax_list) {$tax_levels{$_}=1}

if (not $output_file) {$output_file = $assembly_file . ".txt"};

print "Usage: gc_cov_annotate.pl --blasttaxid CONTIGTAXIDFILE --assembly ASSEMBLYFASTAFILE [--taxdump TAXDUMPDIR] [--bam BAMFILE...] [--cov COVFILES...] [--taxlist species...]\n" .
    "--taxdump is '.' by default, i.e. the files nodes.dmp and names.dmp from the NCBI taxonomy database ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz are expected to be in the current directory\n" . 
    "--taxlist: default is: species order phylum superkingdom, but you can add any other NCBI taxlevel such as class family suborder etc\n" unless
    (-r $blasttaxid_file or $blasttaxid_file eq "-") and -r "$taxdump_dir/nodes.dmp" and -r "$taxdump_dir/names.dmp" and 
    (-r $assembly_file or $assembly_file eq "-");

die "Please specify --blasttaxid\n" unless (-r $blasttaxid_file or $blasttaxid_file eq "-");
die "Please specify --taxdump\n" unless -r "$taxdump_dir/nodes.dmp" and -r "$taxdump_dir/names.dmp";
die "Please specify --assembly\n" unless (-r $assembly_file or $assembly_file eq "-");

#-----------------------------------------------------------------
# Get taxon annotation info from contig-taxid file
#-----------------------------------------------------------------

my (%taxid_has_parent, %taxid_has_taxlevel, %taxid_has_name);
print STDERR scalar localtime() . " - Loading $taxdump_dir/nodes.dmp and $taxdump_dir/names.dmp into memory ...\n";
&load_nodes_names ("$taxdump_dir/nodes.dmp","$taxdump_dir/names.dmp");
print STDERR scalar localtime() . " - Loading $taxdump_dir/nodes.dmp and $taxdump_dir/names.dmp into memory ... DONE\n";

my $blasttaxid_fh = &read_fh($blasttaxid_file);
my %contig_taxinfo;
while (<$blasttaxid_fh>) {
    die "Contig-taxid file $blasttaxid_file does not seem to have two cols with the seqid in the first col and taxid in the second col" unless 
        /^(\S+)\t(\d+)/;
    $contig_taxinfo{$1} = &taxonomy_report($2);
}
close $blasttaxid_fh;

#-----------------------------------------------------------------
# Calculate GC, len, for assembly file, get coverage from bamfiles
#-----------------------------------------------------------------

# load assembly file into memory:
print STDERR scalar localtime() . " - Loading assembly fasta file $assembly_file into memory ...\n";
my $fastahash = &fastafile2hash ($assembly_file);
print STDERR scalar localtime() . " - Loading assembly fasta file $assembly_file into memory ... DONE\n";

# calculate coverage for each bam file
for my $bam_file (@bam_files) {
    my $reads_processed = 0;
    my ($refstart, $refend, $readid, @F); # declaring here to avoid having to declare in each instance of loop
    open SAM, "samtools view $bam_file |" or die $!; 
    print STDERR scalar localtime() . " - Reading $bam_file ...\n";
    while (<SAM>) {
        # ignore comment lines/ SAM headers
        next if /^\s*$/ or /^@/ or /^#/;
        chomp;
        @F=split/\t/;
        if (($F[1] & 4) != 4) {
            $refstart =  $F[3];
            $refend   =  $F[3] - 1;
            while ($F[5] =~ /(\d+)[MIDNP]/g) { $refend += $1 } # calculate coverage on reference sequence from CIGAR
            die "ContigID $F[2] in bamfile $bam_file, but not in assembly file $assembly_file\n" if not exists $$fastahash{$F[2]};
            $$fastahash{$F[2]}{$bam_file} += ($refend - $refstart + 1)/$$fastahash{$F[2]}{len};
        }
        $reads_processed++;
        print STDERR "Processed $reads_processed reads by " . localtime() . "...\n" if $reads_processed % 1000000 == 0;
    }
    print STDERR scalar localtime() . " - Reading $bam_file ... DONE\n";
}

# calculate coverage for each cov file (two col file with seqid in first col, average read coverage depth in second col)
for my $cov_file (@cov_files) {
    open COV, $cov_file or die $!;
    print STDERR scalar localtime() . " - Reading $cov_file ...\n";
    while (<COV>) {
        /^(\S+)\s+(\S+)/ and $$fastahash{$1}{$cov_file} = $2;
    }
    print STDERR scalar localtime() . " - Reading $cov_file ... DONE\n";
}

# print len gc cov (one for each bam, and each cov) and taxon annotation info in one large file (that can be used for plotting by R)
print STDERR scalar localtime() . " - Making len gc cov taxon annotation data file $assembly_file.txt for R ...\n";
open  LENCOVGC, ">$output_file" or die $!;

# header row:
print LENCOVGC "seqid\tlen\tgc";
foreach (@bam_files) {print LENCOVGC "\tcov_$_"};
foreach (@cov_files) {print LENCOVGC "\tcov_$_"};
foreach (@tax_list)  {print LENCOVGC "\ttaxlevel_$_"};
print LENCOVGC "\n";

# for each contig/seqid:
my ($seqid, $length, $gccount, $nonatgc, $totalcov, $cov); # declared outside loop for efficiency
for $seqid (keys %{$fastahash}) { 
    $length   = $$fastahash{$seqid}{len};
    $gccount  = $$fastahash{$seqid}{gc};
    $nonatgc  = $$fastahash{$seqid}{nonatgc};
    print LENCOVGC "$seqid\t$length\t". ($gccount/($length-$nonatgc));
    for my $bam_file (@bam_files) {
        $cov = (exists($$fastahash{$seqid}{$bam_file}) ? $$fastahash{$seqid}{$bam_file} : 0);
        print LENCOVGC "\t" . $cov;
    }
    for my $cov_file (@cov_files) {
        $cov = (exists($$fastahash{$seqid}{$cov_file}) ? $$fastahash{$seqid}{$cov_file} : 0);
        print LENCOVGC "\t" . $cov;
    }
    for my $tax_level (@tax_list) {
        print LENCOVGC "\t" . (exists(${$contig_taxinfo{$seqid}}{$tax_level}) ? ${$contig_taxinfo{$seqid}}{$tax_level} : "Not annotated"); 
    }
    print LENCOVGC "\n";
}
print STDERR scalar localtime() . " - Making len gc cov taxon annotation data file $assembly_file.txt for R ... DONE\n";

############################################################

sub taxonomy_report {
    my $hit_taxid = shift @_;
    my @parents = &get_parents($hit_taxid);
    # convert @parents to tax names:
    my %taxinfo;
    # my $taxonomy_report_string = "";
    for my $parent (@parents) {
        if (exists $taxid_has_taxlevel{$parent} and exists $tax_levels{$taxid_has_taxlevel{$parent}}) {
            $taxinfo{$taxid_has_taxlevel{$parent}} = $taxid_has_name{$parent};
        }
    }
    return \%taxinfo;
    # for my $tax_level (keys %hit_counts) {
        # for my $tax_name (keys %{$hit_counts{$tax_level}}) {
            # $taxonomy_report_string .= "$tax_level\t$tax_name\t";
        # }
    # }
    # return $taxonomy_report_string . "\n";
}

############################################################

sub get_parents {
    my @all = @_;
    my $current_id = $all[0];
    if (exists $taxid_has_parent{$current_id} and $current_id ne $taxid_has_parent{$current_id}) {
        unshift @all, $taxid_has_parent{$current_id};
        @all = &get_parents(@all);
    }
    return @all;
}

############################################################

sub load_nodes_names {
    my $fh;
    my $nodesfile = shift @_;
    my $namesfile = shift @_;
    $fh = &read_fh($nodesfile);
    while (my $line = <$fh>) {
        # line in nodes.dmp should match the regexp below.
        # Change the regexp if NCBI changes their file format
        next if $line !~ /^(\d+)\s*\|\s*(\d+)\s*\|\s*(.+?)\s*\|/;
        $taxid_has_parent{$1} = $2;
        $taxid_has_taxlevel{$1} = $3;
    }
    close $fh;
    
    $fh = &read_fh($namesfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^(\d+)\s*\|\s*(.+?)\s*\|.+scientific name/;
        $taxid_has_name{$1} = $2;
    }
}

############################################################

sub fastafile2hash {
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $seqid;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $seqid = $1;
            $sequences{$seqid}{desc} = $2;
        }
        else {
            chomp $line;
            $sequences{$seqid}{seq}     .= $line;
            $sequences{$seqid}{len}     += length $line;
            $sequences{$seqid}{gc}      += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences{$seqid}{nonatgc} += ($line =~ tr/N/N/);
        }
    }
    close $fh;
    return \%sequences;
}

############################################################

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

############################################################
