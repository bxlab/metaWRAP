#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($file, $remove, $delim) = ("", 0, "\t");
GetOptions (
  "file:s" => \$file,
  "v|remove" => \$remove,
  "delim:s" => \$delim,
);

die "Usage: fgrep.pl [-v] [-d <delimiter] -f <filewithpatterns> <file to be searched>\n" unless $file;

open FILE, "<$file" or die "Could not open file with patterns\n";
my %hash;
while (<FILE>)
{
	chomp;
	$hash{$_}++
}
while (my $line = <>)
{
	chomp $line;
	my $found = 0;
	foreach (split(/$delim/,$line))
	{
		 if (exists $hash{$_})
		 {
		 	$found = 1;
		 	last
		 }
	}
	print "$line\n" if $found and not $remove;
	print "$line\n" if not $found and $remove;
}
