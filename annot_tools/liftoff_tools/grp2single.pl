#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  for my $tb (@ta[2..$#ta]) {
    print "$tb\t$ta[0]\n";
  }
}
