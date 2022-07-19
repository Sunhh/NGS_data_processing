#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 in.gff3 > simple_IDs.gff3\n";

while (<>) {
  m!^\s*(#|$)! and do { print; next; };
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[2] =~ m!^gene$!i) {
    $ta[8] =~ m!(?:^|\s|;)ID=([^\s;]+)! or die "$_\n";
    $ta[8] = "ID=$1";
  } elsif ($ta[2] =~ m!^mRNA$!i) {
    my ($mid, $gid);
    $ta[8] =~ m!(?:^|\s|;)ID=([^\s;]+)! or die "$_\n";
    $mid = $1;
    $ta[8] =~ m!(?:^|\s|;)Parent=([^\s;]+)!i and $gid = $1;
    if (defined $gid) {
      $ta[8] = "ID=$mid;Parent=$gid";
    } else {
      $ta[8] = "ID=$mid";
    }
  } elsif ($ta[2] =~ m!^CDS|exon$!i) {
    $ta[8] =~ m!(?:^|\s|;)Parent=([^\s;]+)!i or die "$_\n";
    $ta[8] = "Parent=$1";
  } else {
    die "[Err] Unknown feature [$ta[2]]\n";
  }
  print join("\t", @ta)."\n";
}
#21QDX551_Chr02	EVM	gene	77439	81626	.	+	.	ID=21QDX551C02G000010;Name=EVM%20prediction%2021QDX551_Chr02.1
#21QDX551_Chr02	EVM	mRNA	77439	81626	.	+	.	ID=21QDX551C02G000010.1;Parent=21QDX551C02G000010;Name=EVM%20prediction%2021QDX551_Chr02.1
#21QDX551_Chr02	EVM	exon	77439	77768	.	+	.	ID=evm.model.21QDX551_Chr02.1.exon1;Parent=21QDX551C02G000010.1

