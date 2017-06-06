#!/usr/bin/env perl

use strict;
use warnings;

our $USAGE = "$0 contigOverlaps.txt firstBinsContigNames.fofn secondBinsContigNames.fofn [...]\n";

die $USAGE unless scalar(@ARGV >= 3);

my $contigsOverlapFile = shift;

open(my $overlapfh, "<", $contigsOverlapFile) || die;
my $overlapHeader = <$overlapfh>;

our @llhh_assemAassemBbinAbinBCounts;
our @llh_assemBinReadCounts;
our @lh_assemContigBins;

foreach my $binfofn (@ARGV) {
    open(my $fhfofn, "<", $binfofn) || die "Could not open $binfofn! $!";
    my %h_contigBins;
    push @lh_assemContigBins, \%h_contigBins;
    while (my $binFile = <$fhfofn>) {
        chomp($binFile);
        open(my $fh, "<", $binFile) || die "Could not open $binFile! $!";;
        my $binName = $binFile;
        $binName =~ s#.*/##;

        while (my $contigName = <$fh>) {
            chomp($contigName);
            $contigName =~ s/[ \t].*//;

            # hack to fix reference name mangling
            $contigName =~ s/[\[]/_/;
            $contigName =~ s/[\]]//;

            $h_contigBins{$contigName} = $binName;
#print STDERR "$contigName -> $binName\n";
        }
        close($fh) || die;
    }
    close($fhfofn) || die
    
}

while (<$overlapfh>) {
    my ($assemblyi, $assemblyj, $contigi, $contigj, $overlap, $totali, $percent) = split(/\t/);
    $contigi =~ s/\s.*//;
    $contigj =~ s/\s.*//;
    
    if (exists $lh_assemContigBins[$assemblyi]{$contigi}) {
      my $binAName = $lh_assemContigBins[$assemblyi]{$contigi};
      my $binBName = "*";
      if ($binBName ne $contigj) {
          $binBName = $lh_assemContigBins[$assemblyj]{$contigj};
          if (not defined $binBName) {
#warn "There is no bin for $contigj in $ARGV[$assemblyj]\n";
              next;
          }
      }
    
      $llh_assemBinReadCounts[$assemblyi][$assemblyj]{$binAName} += $overlap;
      $llhh_assemAassemBbinAbinBCounts[$assemblyi][$assemblyj]{$binAName}{$binBName} += $overlap;
    }
    
}

printf("assemblyA\tassemblyB\tbinA\tbinB\toverlap\ttotalA\t%%\n");
for(my $assemblyi = 0; $assemblyi < scalar(@ARGV); $assemblyi++) {
    for(my $assemblyj = 0; $assemblyj < scalar(@ARGV); $assemblyj++) {
        my $rh_binReadCounts = $llh_assemBinReadCounts[$assemblyi][$assemblyj];

        if ($assemblyi == $assemblyj) { next; }
        my $rh_BinABinBCounts = $llhh_assemAassemBbinAbinBCounts[$assemblyi][$assemblyj];
        foreach my $binAName (sort keys %{$rh_binReadCounts}) {
            my $total = $rh_binReadCounts->{$binAName};
            while(my($binBName, $overlap) = each %{$llhh_assemAassemBbinAbinBCounts[$assemblyi][$assemblyj]{$binAName}}) {
                my $frac = $total > 0 ? $overlap/$total : 0;
                if ($frac > 0.001) {
                  printf("%d\t%d\t%s\t%s\t%d\t%d\t%0.2f\n", $assemblyi, $assemblyj, $binAName, $binBName, $overlap, $total, 100*$frac); 
                }
            }
        }
    }
}
        
