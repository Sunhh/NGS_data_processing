#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

!@ARGV and die "perl $0 s0.anchors.gff3 toAdd_CDS.gff3 out_prefix\n";

my $f1 = shift;
my $f2 = shift;
my $opre = shift;

open F1,'<',"$f1" or die;
my %loc1;
while (<F1>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  $ta[2] eq 'gene' or next;
  push(@{$loc1{$ta[0]}}, [@ta[3,4]]);
}
close F1;
for (keys %loc1) { @{$loc1{$_}} = sort {$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$loc1{$_}}; }

open F2,'-|',"perl ~/tools/github/NGS_data_processing/temp/deal_gff3.pl -inGff $f2 -getJnLoc " or die;
&fileSunhh::write2file("${opre}.list","", '>');
while (<F2>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'mrnaID' and next;
  my $is_ovl = 0;
  for my $a1 (@{$loc1{$ta[2]}}) {
    $a1->[1] < $ta[6] and next;
    $a1->[0] > $ta[7] and last;
    $is_ovl = 1;
    last;
  }
  $is_ovl == 1 and next;
  &fileSunhh::write2file("${opre}.list", "$ta[0]\n",'>>');
}
close F2;
&runCmd("perl ~/tools/github/NGS_data_processing/temp/deal_gff3.pl -inGff $f2 -gffret ${opre}.list -idType mRNA > ${opre}.gff3");

