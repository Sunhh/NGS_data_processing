#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0  in_go.obo.tab  in_GO.annot  out_file\n";

my $ifnTab = shift;
my $ifnAnn = shift;
my $ofnGOen = shift;

my (%go2ancestors, %go2names, %goAlt2Repre);

open F1,'<',"$ifnTab" or die;
while (<F1>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'GO_ID' and next;
  $go2ancestors{$ta[0]} = [split(/;/, $ta[5])];
  $go2names{$ta[0]}{'root'} = $ta[2];
  $go2names{$ta[0]}{'name'} = $ta[3];
  for my $id1 (split(/;/, $ta[1])) {
    $id1 eq 'NA' and next;
    $id1 eq '' and next;
    $goAlt2Repre{$id1} = $ta[0];
  }
}
close F1;

open F2,'<',"$ifnAnn" or die;
open O1,'>',"$ofnGOen" or die; # geneID \t OG_ID \t OG_root \t OG_name \n
my %hasOut;
my %notFound;
while (<F2>) {
  chomp;
  my ($geneID, $goID)=split(/\t/, $_);
  $goID //= "";
  if ($goID eq "") {
    print O1 join("\t", $geneID, "", "", "")."\n";
    $hasOut{$geneID}{""} = 1;
  }
  defined $goAlt2Repre{$goID} and $goID = $goAlt2Repre{$goID};
  if (!(defined $go2ancestors{$goID})) {
    $notFound{$geneID}{$goID} = 1; next;
  }
  for my $id1 (@{$go2ancestors{$goID}}) {
    defined $hasOut{$geneID}{$id1} and next;
    print O1 join("\t", $geneID, $id1, $go2names{$id1}{'root'}, $go2names{$id1}{'name'})."\n";
    $hasOut{$geneID}{$id1} = 1;
  }
}

for my $geneID (sort keys %notFound) {
  defined $hasOut{$geneID} and next;
  print O1 join("\t", $geneID, "", "", "")."\n";
  # for my $id1 (sort keys %{$notFound{$geneID}}) {
  #   print O1 join("\t", $geneID, $id1, "", "")."\n";
  # }
}
close O1;
close F2;

