#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0  synOG.grp  label_DEGs.tsv > label_DEG_groups.tsv\n";

my $f1grp = shift;
my $f2deg = shift; # 'N|U|D'

my %gene2grp;
open F1,'<',"$f1grp" or die;
while (<F1>) {
  chomp;
  my @ta=split(/\s|,/, $_);
  for (my $i=1; $i<@ta; $i++) {
    $ta[$i] =~ m!^\s*$! and next;
    $gene2grp{$ta[$i]} = $ta[0];
  }
}
close F1;

my (@grp_lines, %grp_lineIdx, @grp_names);
open F2,'<',"$f2deg" or die;
while (<F2>) {
  chomp;
  my @ta=split(/\t/, $_);
  if ($. == 1) {
    print STDOUT join("\t", @ta)."\n";
    next;
  }
  my ($geneID) = shift(@ta);
  $gene2grp{$ta[0]} //= "$geneID";
  my $grpID = $gene2grp{$geneID};
  $grp_lineIdx{$grpID} //= scalar(@grp_lines);
  $grp_names[$grp_lineIdx{$grpID}] //= $grpID;
  push(@{$grp_lines[$grp_lineIdx{$grpID}]}, [@ta]);
}
close F2;


for (my $i=0; $i<@grp_names; $i++) {
  my @o = ($grp_names[$i]);
  for (my $j=0; $j<@{$grp_lines[$i][0]}; $j++) {
    my %label_h;
    for (my $m=0; $m<@{$grp_lines[$i]}; $m++) {
      $label_h{$grp_lines[$i][$m][$j]} ++;
    }
    delete($label_h{'N'});
    my $grp_label = (scalar(keys %label_h) > 0) ? join(";", sort {$label_h{$b} <=> $label_h{$a} || $b cmp $a} keys %label_h) : 'N' ;
    push(@o, $grp_label);
  }
  print STDOUT join("\t", @o)."\n";
}

