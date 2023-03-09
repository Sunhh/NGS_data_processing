#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 num_chr min_monomer_len whole/after_clusterin_summary_whole.txt.cntChr > whole/after_clusterin_summary_whole.txt.cntChr.slct\n";

my $nChr = shift;
my $min_monomer_len = shift;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq "CL_ID" and next;
  $ta[4] == $nChr or next;
  $ta[2] >= $min_monomer_len or next;
  my @a1 = split(/;/, $ta[5]);
  my @a2 = split(/;/, $ta[6]);
  for (my $i=0; $i<@a1; $i++) {
    my @a3 = split(/\<\*\>/, $a2[$i]);
    print STDOUT join("\t", $a2[$i], "$a1[$i]__$ta[0]", $a1[$i], $ta[0], $a3[2], $a3[3])."\n";
  }
}

