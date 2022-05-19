#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;

!@ARGV and die "perl $0  CLpan  CLpan.trim2CDS.gff3.JnLoc > CLpan.trim2CDS.blk\n";

my $tag = shift;
my $fn_JnLoc = shift;

open F,'<', $fn_JnLoc or die;
while (<F>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'mrnaID' and next;
  print STDOUT join("\t", "$tag:$ta[0]", @ta[2,5,9])."\n";
}
close F;

