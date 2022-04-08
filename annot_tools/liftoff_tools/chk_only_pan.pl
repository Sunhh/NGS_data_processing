#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 comb.grp2.novl_loc.wiRepre.fmt > comb.grp2.novl_loc.wiRepre.fmt.ifOnlyPan\n";

while (<>) {
  chomp;
  $. == 1 and do { print "$_\n"; next; };
  my @ta=split(/\t/, $_);
  my @ifpan;
  for (my $i=7; $i<@ta; $i++) {
    for my $tb (split(/;/, $ta[$i])) {
      $tb =~ m!^\S+:\d+\-\d+:[+-]$! and next;
      $tb =~ m!^C.Pan! and next;
      $ifpan[$i-7] ++;
    }
  }
  for (my $i=0; $i<4; $i++) {
    $ifpan[$i] //= 0;
  }
  print join("\t", @ta[0..6], @ifpan)."\n";
}
