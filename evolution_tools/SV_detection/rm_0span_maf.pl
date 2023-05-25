#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 align2w38_anc_fix_jn.maf > fixed.maf\n";

while (<>) {
  chomp;
  if (m!^\s*#!) {
    print STDOUT "$_\n";
    next;
  }
  my $l2=<>; chomp($l2);
  my $l3=<>; chomp($l3);
  my $l4=<>; chomp($l4);
  m!^a\s! or die "bad2:$_\n";
  $l2 =~ m!^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s! or die "bad3:$l2\n";
  $3 == 0 and next;
  $l3 =~ m!^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s! or die "bad4:$l3\n";
  $3 == 0 and next;
  print "$_\n$l2\n$l3\n$l4\n";
}

