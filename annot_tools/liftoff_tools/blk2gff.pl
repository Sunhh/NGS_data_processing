#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 output/CX.cds.blk.CL > output/CX.cds.blk.CL.gff3\n";

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my $i=0;
  my (@cds, $cS, $cE);
  for my $tb (split(/;/, $ta[3])) {
    $tb =~ m!^(\d+)\,(\d+)$! or die "$tb\n";
    push(@cds, [$ta[1], "blk", "CDS", $1, $2, ".", $ta[2], ".", "Parent=$ta[0]"]);
    $cS //= $1; $cS > $1 and $cS = $1;
    $cE //= $1; $cE < $2 and $cS = $2;
  }
  print join("\t", $ta[1], "blk", "gene", $cS, $cE, ".", $ta[2], ".", "ID=$ta[0]-G")."\n";
  print join("\t", $ta[1], "blk", "mRNA", $cS, $cE, ".", $ta[2], ".", "ID=$ta[0];Parent=$ta[0]-G")."\n";
  for my $a1 (@cds) {
    print join("\t", @$a1)."\n";
  }
}

