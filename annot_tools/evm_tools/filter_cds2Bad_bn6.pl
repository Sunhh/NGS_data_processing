#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 evmMerged.c.fa.toBad.bn6 > evmMerged.c.fa.toBad.bn6.toRM\n";

my $min_ident = 97;
my $min_percCov  = 95; # 

my %h;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  defined $h{$ta[0]} and next;
  $ta[2] >= $min_ident or next;
  $ta[7]-$ta[6]+1 >= $ta[12] * $min_percCov / 100 or next;
  $h{$ta[0]} = 1;
  print STDOUT "$_\n";
}

