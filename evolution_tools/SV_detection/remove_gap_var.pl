#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 flank_len_100 C31.ndf_ref_snps.gff > C31.ndf_ref_snps.rmGap.gff\n";

my $flank_len = shift;

my (@lines, %gaps);
while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[8] =~ m!Name=[^;\s]*(gap|ATGCN)!) {
    push(@{$gaps{$ta[0]}}, [@ta[3,4]]);
  } else {
    push(@lines, [@ta]);
  }
}

for my $cid (keys %gaps) {
  @{$gaps{$cid}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$gaps{$cid}};
}
for my $l1 (@lines) {
  my $close2gap = 0;
  for my $seR (@{$gaps{$l1->[0]}}) {
    $seR->[1] < $l1->[3]-$flank_len and next;
    $seR->[0] > $l1->[4]+$flank_len and last;
    $close2gap = 1;
    last;
  }
  $close2gap == 1 and next;
  print STDOUT join("\t", @$l1)."\n";
}

