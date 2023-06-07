#!/usr/bin/perl
# [3/9/2023] Only applicable after running 'bcftools norm -m-both | bcftools norm -d none --fasta-ref xx.fa | bcftools sort ';
use strict;
use warnings;

-t and !@ARGV and die "perl $0 normSrt_dup.vcf > dedup.vcf\n";

my $prevLine = '';
my ($prevID, $prevPos, $prevLen) = ("", -1, -1);
while (<>) {
  m!^\s*(#|$)! and do {print; next;};
  chomp;
  my @ta=split(/\t/, $_);
  my $currLen = abs(length($ta[3])-length($ta[4]));
  $ta[4] eq '<INV>' and $currLen = 1e9;
  if ($ta[0] eq $prevID and $ta[1] == $prevPos) {
    # dup
    if ($currLen > $prevLen) {
      $prevLine = $_;
      $prevLen = $currLen;
    }
  } else {
    $prevLine ne '' and print STDOUT "$prevLine\n";
    $prevLine = $_;
    $prevID = $ta[0];
    $prevPos = $ta[1];
    $prevLen = $currLen;
  }
}
if ($prevLine ne '') {
  print STDOUT "$prevLine\n";
  $prevLine = '';
}

