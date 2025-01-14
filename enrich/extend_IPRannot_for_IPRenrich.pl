#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 pub_IPR-entry.list  in_gene2IPR.combined  out_file\n";


my $ifnIPRentry = shift;
my $ifhGene2IPR = shift;
my $ofnEnrich   = shift;

my (%ipr2name);

open F1,'<',"$ifnIPRentry" or die;
while (<F1>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'ENTRY_AC' and next;
  $ipr2name{$ta[0]} = $ta[2];
}
close F1;

open F2,'<',"$ifhGene2IPR" or die;
open O1,'>',"$ofnEnrich" or die;
my (%hasOut, %notFound);
while (<F2>) {
  my @ta=split(/\t/, $_);
  chomp(@ta);
  if ($ta[1] eq "" or $ta[1] =~ m!^na$!i) {
    print O1 join("\t", $ta[0], "", "")."\n";
    $hasOut{$ta[0]}{""} = 1;
    next;
  }
  for my $id1 (split(/;/, $ta[1])) {
    if (defined $ipr2name{$id1}) {
      defined $hasOut{$ta[0]}{$id1} and next;
      print O1 join("\t", $ta[0], $id1, $ipr2name{$id1})."\n";
      $hasOut{$ta[0]}{$id1} = 1;
    } else {
      $notFound{$ta[0]}{$id1} = 1;
    }
  }
}
for my $gene1 (sort keys %notFound) {
  defined $hasOut{$gene1} and next;
  print O1 join("\t", $gene1, "", "")."\n";
}
close O1;
close F2;

