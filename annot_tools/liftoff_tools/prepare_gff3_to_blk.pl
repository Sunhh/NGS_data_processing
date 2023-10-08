#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0  CLpan  CLpan.trim2CDS.gff3.JnLoc > CLpan.trim2CDS.blk\n";

my $tag = shift;
my $fn_JnLoc = shift;

$tag ne '' and $tag .= ':';

open F,'<', $fn_JnLoc or die;
while (<F>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'mrnaID' and next;
  my @se = split(/;/, $ta[9]);
  my $cdsL = 0;
  for my $tb (@se) {
    $tb =~ m!^(\d+)\,(\d+)$! or die "[Err] bad blk format |$tb|\n";
    $cdsL += ($2-$1+1);
  }
  print STDOUT join("\t", "$tag$ta[0]", @ta[2,5,9], '.', $cdsL)."\n";
}
close F;

