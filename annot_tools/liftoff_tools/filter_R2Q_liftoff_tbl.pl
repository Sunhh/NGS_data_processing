#!/usr/bin/perl
use strict;
use warnings;

my @a1;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'trans_R_ID' and do { print "$_\n"; next; };
  # (1) Require coverage > 0.9;
  $ta[1] >= 0.9 or next;
  # (2) Require at least 100 bp or 50% of shorter CDS overlapping. or no overlapping to any gene.
  $ta[9] eq '.' or $ta[3] >= 100 or $ta[3] >= 0.5 * $ta[7] or $ta[3] >= 0.5 * $ta[10] or next;
  print "$_\n";
}

