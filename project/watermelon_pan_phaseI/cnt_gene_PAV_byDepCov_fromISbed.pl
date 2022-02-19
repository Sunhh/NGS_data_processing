#!/usr/bin/perl
# Count the PAV of genes in a bed file based on the read depth and coverage information from bedtools output file.
# CMD: bedtools intersect -a ../../06_gene_prediction/final/CLpan.gff3.CDS.bed -b WM78_XMF.CLpan.bam.gcov -wao
use strict;
use warnings;

!@ARGV and die "perl $0 2_rdDep 90_coverage WM9_CalhounGray.CLpan.bam.gcovH.intersect > PAV.list\n";

my $exe_bedtools = 'bedtools';

my $rdD    = shift; # minimum read depth required as a present segment.
my $covP   = shift; # minimum coverage proportion required as a present gene.
my $bedIS  = shift;

my (%cov, %info, @geneList);
open F,'<',"$bedIS" or die;
print STDOUT join("\t", qw/geneID PAV CDS_len CDS_cov CDS_covP/)."\n";
while (<F>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[3] =~ m!^(\S+)_HS(\d+)$! or die "[Err] Bad format in [$bedIS]\n";
  my ($geneID, $partN) = ($1, $2);
  defined $info{$geneID} or push(@geneList, $geneID);
  my $tk = join("\t", @ta[0..3]);
  unless (defined $info{$geneID}{'seg'}{$tk}) {
    $info{$geneID}{'seg'}{$tk} = 1;
    $info{$geneID}{'len'} += $ta[2] - $ta[1];
  }
  $ta[7] eq '.' and next;
  $ta[7] >= $rdD or next;
  $cov{$geneID} += $ta[8];
}
close F;

for my $geneID (@geneList) {
  my $tag = 'A';
  $cov{$geneID} //= 0;
  if ($cov{$geneID} >= $info{$geneID}{'len'} * $covP/100) {
    $tag = "P";
  }
  print STDOUT join("\t", $geneID, $tag, $info{$geneID}{'len'}, $cov{$geneID}, sprintf("%.2f", 100*$cov{$geneID}/$info{$geneID}{'len'}))."\n";
}

# Cla97Chr01      2023    2305    Cla97C01G000010.1_HS1   Cla97Chr01      2045    2050    5       5
# Cla97Chr01      2023    2305    Cla97C01G000010.1_HS1   Cla97Chr01      2050    2065    4       15

