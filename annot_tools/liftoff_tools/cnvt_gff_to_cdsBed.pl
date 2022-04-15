#!/usr/bin/perl
use strict;
use warnings;

my $htxt = <<HH;

perl $0 out_prefix in.gff3

# Result files are:
  out_prefix.p.CDS.bed
  out_prefix.m.CDS.bed

HH

!@ARGV and die "$htxt";
# perl tools/cnvt_gff_to_cdsBed.pl $_.gff3 > $_.CDS.bed

my $opref = shift;

open OP,'>',"$opref.p.CDS.bed" or die;
open OM,'>',"$opref.m.CDS.bed" or die;
my %h;
while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  $ta[2] =~ m!^CDS$!i or next;
  $ta[8] =~ m!Parent=([^\s;]+)! or die "$_\n";
  $h{$1} ++;
  if      ($ta[6] eq '+') {
    print OP join("\t", $ta[0], $ta[3]-1, $ta[4], "${1}_HS$h{$1}", $1)."\n";
  } elsif ($ta[6] eq '-') {
    print OM join("\t", $ta[0], $ta[3]-1, $ta[4], "${1}_HS$h{$1}", $1)."\n";
  } else {
    die "[Err] Bad str line: $_\n";
  }
}

