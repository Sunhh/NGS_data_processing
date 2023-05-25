#!/usr/bin/perl
use strict;
use warnings;


while (<>) {
  chomp;
  m!^\s*(#|$)! and next;
  my @ta=split(/\t/, $_);
  if ($ta[0] eq 'refChr') {
    print join("\t", qw/refSpan qrySpan/, @ta)."\n";
    next;
  }
  my $lenR = $ta[2]-$ta[1]+1;
  my $lenQ = $ta[5]-$ta[4]+1;
  print join("\t", $lenR, $lenQ, @ta)."\n";
}

