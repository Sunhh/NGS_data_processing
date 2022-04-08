#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 comb.grp2 > comb.grp2.novl_loc\n";

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my (@gID, @loci);
  for my $tb (@ta[2..$#ta]) {
    if ($tb =~ m!^(\S+):(\d+)\-(\d+):([+-])$!) {
      push(@loci, [$1, $2, $3, $4, $3-$2+1]); # chrID, start, end, str, span_length
    } else {
      push(@gID, $tb);
    }
  }
  my @new_loci;
  @loci = sort { $b->[4] <=> $a->[4] } @loci;
  for my $l1 (@loci) {
    my $is_ovl = 0;
    for my $l2 (@new_loci) {
      $l1->[0] eq $l2->[0] or next;
      $l1->[3] eq $l2->[3] or next;
      $l1->[1] > $l2->[2] and next;
      $l1->[2] < $l2->[1] and next;
      $is_ovl = 1;
      last;
    }
    $is_ovl == 0 and push(@new_loci, [@$l1]);
  }
  print STDOUT join("\t", $ta[0], scalar(@gID)+scalar(@new_loci), @gID, (map { "$_->[0]:$_->[1]-$_->[2]:$_->[3]" } @new_loci))."\n";
}

