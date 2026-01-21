#!/usr/bin/perl
use strict;
use warnings;

my $help_txt = <<HH;
######################################################
# perl $0 list-columns mat.grp_melt > mat.grp
# 
#   Format of grp_melt: OGID \\t cellValue \\t columnID
#
# # The missing value will be empty.
######################################################

HH

!@ARGV and die "$help_txt";

my $fnCol = shift;
my $missVal = "";

open F1,'<',"$fnCol" or die;
my (@cols);
while (<F1>) { chomp; push(@cols, (split)[0]); }
close F1;

my %ogRank;
my %cellFill;
while (<>) {
  chomp;
  my ($ogID, $cellVal, $colID) = split(/\t/, $_);
  $ogID eq "OGID" and next;
  $ogRank{$ogID} //= $.;
  push(@{$cellFill{$ogID}{$colID}}, $cellVal);
}
print STDOUT join("\t", "OGID", @cols)."\n";
for my $ogID (sort {$ogRank{$a} <=> $ogRank{$b}} keys %cellFill) {
  my @o1=($ogID);
  for my $colID (@cols) {
    $cellFill{$ogID}{$colID} //= [];
    push(@o1, join(", ", @{$cellFill{$ogID}{$colID}}));
  }
  print STDOUT join("\t", @o1)."\n";
}

