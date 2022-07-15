#!/usr/bin/perl
use strict; 
use warnings;

my @rec; 
my ($geneID, $mrnaID); 
my %usedID; 
my $has_exon = 0; 

my $cnt = 0; 
$cnt++; 

my @badID = qw(
CsGy2G027260
CsGy2G022560
CsGy2G016410
CsGy2G019230
CsGy5G030020
);

LINE: 
while (<>) { 
  m!^\s*(#|$)! and next; 
  chomp; 
  my @ta=split(/\t/, $_); 
  for my $tormID (@badID) {
    $ta[8] =~ m!$tormID! and next LINE; 
  }
  if ($ta[2] eq 'gene') {
    if (scalar(@rec) > 0) {
      &chkout(); 
    }
    @rec = (); 
    $ta[8] =~ m!^ID=([^\s;]+?)\s*(;|$)! or die "$_\n"; 
    $geneID = $1; 
    $mrnaID = ""; 
    push(@rec, [@ta]); 
    next; 
  }
  if ($ta[2] =~ m!^mRNA$!i) {
    $ta[8] =~ m!^ID=([^\s;]+?)\s*(;|$)! or die "$_\n"; 
    my $cur_mrnaID = $1; 
    $ta[8] =~ m!Parent=([^\s;]+?)\s*(;|$)! or die "mrnaParent: $_\n"; 
    my $cur_geneID = $1; 
    if (scalar(@rec) == 1 and $cur_geneID eq $geneID) {
      $mrnaID = $cur_mrnaID; 
      push(@rec, [@ta]); 
    } else {
      &chkout(); 
      $mrnaID = $cur_mrnaID; 
      $geneID = $cur_geneID; 
      push(@rec, [@ta[0,1], "gene", @ta[3,4,5,6,7], "ID=$cur_geneID"]); 
      push(@rec, [@ta]); 
    }
    next; 
  }

  if ($ta[2] =~ m!^exon$!i) {
    $ta[8] =~ m!Parent=([^\s;]+?)\s*(;|$)! or die "exonParent: $_\n";
    $mrnaID eq $1 or die "diff mrnaID ($1 VS. $mrnaID):$_\n"; 
    $has_exon = 1; 
    push(@rec, [@ta]); 
  }
  if ($ta[2] =~ m!^CDS$!i) {
    $ta[8] =~ m!Parent=([^\s;]+?)\s*(;|$)! or die "cdsParent: $_\n";
    $mrnaID eq $1 or die "diff mrnaID ($mrnaID):$_\n"; 
    push(@rec, [@ta]); 
  }

}

sub chkout {
  for my $la1 (@rec) {
    $cnt ++; 
    if ($la1->[2] =~ m!^(gene|mRNA)$!i) {
      print join("\t", @$la1)."\n"; 
    } elsif ($la1->[2] =~ m!^CDS$!) {
      $has_exon == 0 and print join("\t", @{$la1}[0,1], "exon", @{$la1}[3,4,5,6,7], "ID=hsE$cnt;Parent=$mrnaID\n");
      print join("\t", @{$la1}[0..7], "ID=hsC$cnt;Parent=$mrnaID\n");
    } elsif ($la1->[2] =~ m!^exon$!) {
      print join("\t", @{$la1}[0..7], "ID=hsE$cnt;Parent=$mrnaID\n");
    }
  }
  @rec = (); 
  $has_exon = 0; 
  $geneID = ''; 
  $mrnaID = ''; 
}
