#!/usr/bin/perl
# Randomly choose two alleles from four alleles.
# CMD: bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP --threads 30 data/4x_to_tzn1hap.vcf.gz|perl tools/cnvt_ploidy_4to2.pl| bgzip -@ 30 -c > data/2x.vcf.gz
use strict;
use warnings;

-t and !@ARGV and die "bcftools mpileup -f ref_genome.fa  sample.bam  | bcftools call -c -Ov --ploidy 1 -p 1 | perl $0 > sample-g.vcf\n";

while (<>) {
  m!^#! and do {print; next;};
  chomp;
  my @ta=split(/\t/, $_);
  # $ta[8] =~ m!^GT! or die "Bad format: |$ta[8]| $_\n";
  $ta[8] eq 'GT:PL' or die "Bad format: |$ta[8]| $_\n";
  $ta[8] = "GT";
  for my $b1 (@ta[9..$#ta]) {
    # $b1 =~ s!^([\d.]+)(:|$)!$1/$1$2! or die "|$b1| in line: $_\n";
    $b1 =~ m!^([\d.]+)(:|$)! or die "|$b1| in line: $_\n";
    $b1 = "$1/$1";
  }
  print STDOUT join("\t", @ta)."\n";
}

