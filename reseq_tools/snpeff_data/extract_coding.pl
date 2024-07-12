#!/usr/bin/perl
# 7/11/2024: Each line should have only one gene instead of merged genes.
use strict;
use warnings;

my %t2g=qw(whole coding coding coding splice intron intron intron promoter promoter upstream promoter mapgene remnant downstream downstream intergenic inter coding_low coding);
### coding       : whole + coding;
### intron       : splice + intron;
### promoter     : promoter + upstream;
### remnant      : mapgene; They are 'non-97103 remnant'.
### downstream   : downstream;
###   there should be 'non-97103 completely lost' genes too, which locate in super long insertions. They should be added later from other analysis.

my $cN_gene    = 4;
my $cN_simEff  = 7;
my $cN_dist    = 6;
my $cN_oriEff  = 8;

!@ARGV and die "perl $0 coding/intron/promoter/remnant/downstream allInDel-ann_ud2k.tbl2 > allInDel-ann_ud2k.tbl3\n";

my $want = shift; $want = lc($want);
my %has_line;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[0] eq 'CHROM' or $ta[4] eq '.') {
    print STDOUT join("\t", @ta[0..3], $ta[$cN_gene], 'colEffect', @ta[$cN_dist, $cN_oriEff])."\n";
    next;
  }
  defined $t2g{$ta[$cN_simEff]} or die "[Err] Effect [$ta[$cN_simEff]] not defined!\n";
  my $o;
  if ($t2g{$ta[$cN_simEff]} eq $want) {
    $o=join("\t", @ta[0..3], @ta[$cN_gene, $cN_simEff, $cN_dist, $cN_oriEff])."\n";
  } else {
    $o=join("\t", @ta[0..3], qw/. . . ./)."\n";
  }
  defined $has_line{$o} and next;
  $has_line{$o} = 1;
  print STDOUT $o;
}

