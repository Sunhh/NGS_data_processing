#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 data/slct-synOG.grp > data/slct-synOG.trans2gene.txt\n";

print "transcript_id\tgene_id\n";
while (<>) {
  chomp;
  my @ta=split(/\s|,/, $_);
  $ta[0] eq "OGID" and next;
  my ($grpID) = shift(@ta);
  for my $tb (@ta) {
    $tb =~ m!^\s*$! and next;
    print "$tb\t$grpID\n";
  }
}

