#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 in.gff3 > in.rmGene.gff3\n"; 

while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  scalar(@ta) > 6 or next;
  $ta[2] =~ m!^(mRNA|CDS)$!i or next;
  s!Name=[^\s;]+!!;
  print "$_\n";
}
