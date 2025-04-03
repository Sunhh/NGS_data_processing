#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 combined_DEGs.txt > label_DEGs.tsv\n";

my $minFDR  = 0.05;
my $minlog2 = 1;
my $startLineN = 3;
my $startColN = 1;

my @h1;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my @o = ($ta[0]);
  if ($. == 1) {
    @h1=@ta;
    for (my $i=1;$i<@ta;$i+=4) {
      push(@o, $ta[$i]);
    }
    print STDOUT join("\t", @o)."\n";
    next;
  }
  $. < $startLineN and next;
  for (my $i=$startColN; $i<@ta; $i+=4) {
    my $isDEG = ($ta[$i] eq 'NA' or $ta[$i] >= $minFDR) ? 'N' : ($ta[$i+1] > $minlog2) ? 'U' : ($ta[$i+1] < -$minlog2) ? 'D' : 'N';
    push(@o, $isDEG);
  }
  print STDOUT join("\t", @o)."\n";
}

# GeneID  PI296341_flesh_10D_18D  PI296341_flesh_10D_18D  PI296341_flesh_10D_18D  PI296341_flesh_10D_18D  PI296341_flesh_18D_26D  PI296341_flesh_18D_26D  PI296341_flesh_18D_26D  PI296341_fle>
# GeneID  MeanTPM_Baseline        MeanTPM_Treatment       FDR     Log2FoldChange  MeanTPM_Baseline        MeanTPM_Treatment       FDR     Log2FoldChange  MeanTPM_Baseline        MeanTPM_Trea>
# XG0025C01G000010        1.8614251691445 3.40384358267634        0.257661239477431       0.870756771162025       3.40384358267634        1.37448784610886        0.337933502891573       -1.3>
# XG0025C01G000020        129.792051799986        135.816431021726        0.795783587352077       0.065455988016029       135.816431021726        158.258420206539        0.75472632700123    >
# XG0025C01G000030        125.586952010472        75.9892632170719        0.0913563289826816      -0.724819079435612      75.9892632170719        22.009099453772 0.00352649702928227     -1.7>

