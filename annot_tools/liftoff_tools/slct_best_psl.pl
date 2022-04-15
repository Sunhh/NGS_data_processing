#!/usr/bin/perl
# [3/28/2022] CDS blocks are required for further use.
# [3/29/2022] same strand + Coverage >= 0.5 of either + Blat Identity >= 0.95 + Highest match bases
use strict;
use warnings;

!@ARGV and die "perl $0 in.psl > in.psl.gene_pair\n";

my $min_cov = 0.5;
my $min_blatIdent = 0.95;

my %h1;
while (<>) {
  chomp;
  m!^\d+! or next;
  my @ta=split(/\t/, $_);
  $ta[8] eq '+' or next;
  $ta[12]-$ta[11] >= $min_cov * $ta[10] or $ta[16]-$ta[15] >= $min_cov * $ta[14] or next;
  $ta[0] >= $min_blatIdent * ($ta[0] + $ta[1]) or next;
  $h1{$ta[9]} //= [@ta];
  if      ($h1{$ta[9]} < $ta[0]) {
    $h1{$ta[9]} = [@ta];
  } elsif ($h1{$ta[9]} == $ta[0]) {
    if ( $h1{$ta[9]}[12]-$h1{$ta[9]}[11] > $ta[12]-$ta[11] ) {
      $h1{$ta[9]} = [@ta];
    } elsif ( $h1{$ta[9]}[12]-$h1{$ta[9]}[11] == $ta[12]-$ta[11] ) {
      if ( $h1{$ta[9]}[16]-$h1{$ta[9]}[15] > $ta[16]-$ta[15] ) {
        $h1{$ta[9]} = [@ta];
      }
    }
  }
}
for (keys %h1) {
  my $cov1 = ($h1{$_}[12]-$h1{$_}[11])/$h1{$_}[10];
  my $cov2 = ($h1{$_}[16]-$h1{$_}[15])/$h1{$_}[14];
  print join("\t", $_, $h1{$_}[13], $h1{$_}[0]/($h1{$_}[0] + $h1{$_}[1]), $h1{$_}[0], $cov1, $cov2, $h1{$_}[10], $h1{$_}[14])."\n";
}

