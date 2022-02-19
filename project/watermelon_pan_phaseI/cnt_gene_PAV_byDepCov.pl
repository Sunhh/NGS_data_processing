#!/usr/bin/perl
# Count the PAV of genes in a bed file based on the read depth and coverage information from bedtools output file.
# CMD: bedtools intersect -a ../../06_gene_prediction/final/CLpan.gff3.CDS.bed -b WM78_XMF.CLpan.bam.gcov -wao
use strict;
use warnings;

!@ARGV and die "perl $0 2_rdDep 90_coverage HS_CDS_loc.bed rdCov.bed > PAV.list\n";

my $exe_bedtools = 'bedtools';

my $rdD = shift;
my $covP = shift;
my $bedRef = shift;
my $bedCov = shift;

my (%cov, %info, @geneList);
# Get gene length
open F1,'<', $bedRef or die;
while (<F1>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[3] =~ m!^(\S+)_HS(\d+)$! or die "[Err] Bad format in [$bedRef]\n";
  my ($geneID, $partN) = ($1, $2);
  defined $info{$geneID} or push(@geneList, $geneID);
  defined $info{$geneID}{'seg'}{$ta[3]} and die "repeat!$_\n";
  $info{$geneID}{'len'} += $ta[2] - $ta[1];
}
close F1;

open F,'-|',"$exe_bedtools intersect -a $bedRef -b $bedCov -wao" or die;
print STDOUT join("\t", qw/geneID PAV CDS_len CDS_cov CDS_covP/)."\n";
while (<F>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[7] >= $rdD or next;
  $ta[3] =~ m!^(\S+)_HS(\d+)$! or die "[Err] Bad format in [$bedRef]\n";
  my ($geneID, $partN) = ($1, $2);
  $cov{$geneID} += $ta[8];
}
close F;
for my $geneID (@geneList) {
  my $tag = 'V';
  $cov{$geneID} //= 0;
  if ($cov{$geneID} >= $info{$geneID}{'len'} * $covP/100) {
    $tag = "P";
  }
  print STDOUT join("\t", $geneID, $tag, $info{$geneID}{'len'}, sprintf("%.2f", 100*$cov{$geneID}/$info{$geneID}{'len'}))."\n";
}

# Cla97Chr01      2023    2305    Cla97C01G000010.1_HS1   Cla97Chr01      2045    2050    5       5
# Cla97Chr01      2023    2305    Cla97C01G000010.1_HS1   Cla97Chr01      2050    2065    4       15

