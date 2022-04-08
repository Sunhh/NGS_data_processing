#!/usr/bin/perl
# [3/25/2022] The column number may be different!!!
use strict;
use warnings;

-t and !@ARGV and die "perl $0 bedtools_intersect_wao.out > bedtools_intersect_wao.out.tbl\n";

my $cN_mID_1 = 4;
my $cN_mID_2 = 9;
my $cN_ovlLen = 10;

my %h;
my %trans_cdsLen;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $h{$ta[$cN_mID_1]}{$ta[$cN_mID_2]} += $ta[$cN_ovlLen];
  $trans_cdsLen{$ta[$cN_mID_1]} += $ta[$cN_ovlLen];
}
my @o1;
for my $g1 (keys %h) {
  my @tg2 = keys %{$h{$g1}};
  for my $g2 (@tg2) {
    if ($g2 eq "." and scalar(@tg2) > 1) {
      next;
    }
    push(@o1, [$g1, $g2, $h{$g1}{$g2}]); # [mID_1, mID_2, overlap_size]
  }
}
print STDOUT join("\t", qw/trans_mrnaID tgt_mrnaID trans_cdsLen trans_ovlLen/)."\n";
for my $a1 (sort {$b->[2] <=> $a->[2]} @o1) {
  print STDOUT join("\t", $a1->[0], $a1->[1], $trans_cdsLen{$a1->[0]}, $a1->[2])."\n";
}

