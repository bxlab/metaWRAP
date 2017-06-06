#!/usr/bin/perl

use warnings;
use strict;

our $USAGE = "$0 depth.txt bin1.fa [...]\n\n";
die $USAGE unless scalar(@ARGV) >= 2;

our %hdepths;

open(my $d, "<", $ARGV[0]) || die;
my $header = <$d>;
my $reg = qr/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/;
my $reg2 = qr/^(\S+)/;
while (<$d>) {
  my($contig, $len, $totalDepth) = $_ =~ $reg;
  my($contigName) = $contig =~ $reg2;
  $hdepths{$contigName} = [$len, $totalDepth];
}

close($d) || die;
shift;

printf("bin\ttotalLength\tAvgDepth\tStdDev\n");

my %hbins;
my @lbins;
foreach my $file (@ARGV) {
  open(my $f, "<", $file) || die;
  $reg = qr/^>(\S+)/;
  my @l_contigs;
  $hbins{$file} = \@l_contigs;
  push @lbins, $file;
  while (<$f>) {
    if ($_ =~ $reg) {
      my $contig = $1;
      push @l_contigs, $contig;
    }
  }
  close($f) || die;

  my $totalLen = 0;
  my $totalDepth = 0;

  foreach my $contig (@l_contigs) {
    my $rl = $hdepths{$contig};
    $totalLen += $rl->[0];
    $totalDepth += $rl->[0] * $rl->[1];
  }

  
  my($avgDepth, $stdDev) = (0,0);
  if ($totalLen > 0) {
    $avgDepth = $totalDepth / $totalLen;
    foreach my $contig (@l_contigs) {
      my $rl = $hdepths{$contig};
      my $diff = $avgDepth - $rl->[1];
      $stdDev += $diff * $diff * $rl->[0] / $totalLen;
    }
    $stdDev = sqrt( $stdDev );
  }

  printf("%s\t%d\t%0.2f\t%0.2f\n", $file, $totalLen, $avgDepth, $stdDev);
}

