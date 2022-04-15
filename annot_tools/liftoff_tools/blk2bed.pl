#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 output/CX.cds.blk.CL > output/CX.cds.blk.CL.bed\n";

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my $i=0;
  for my $tb (split(/;/, $ta[3])) {
    $i++;
    $tb =~ m!^(\d+)\,(\d+)$! or die "$tb\n";
    print join("\t", $ta[1], $1-1, $2, "$ta[0]_HS$i")."\n";
  }
}
