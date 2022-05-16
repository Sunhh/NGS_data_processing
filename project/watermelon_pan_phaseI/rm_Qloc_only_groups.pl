#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0  ov1/comb.grp2.novl_loc.syn.grp_tbl > ov1/comb.grp2.novl_loc.syn.grp_tbl.filtered\n";

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my $is_good = 0;
  for my $tb (@ta[2..$#ta]) {
    $tb =~ m!^\S+:\d+\-\d+:[+-]$! and next;
    $is_good = 1;
  }
  $is_good == 1 or next;
  print "$_\n";
}

# cat ov1/comb.grp2.novl_loc.syn.grp_tbl | perl -e 'while (<>) { chomp; m!\t\S+:\d+\-\d+:[+-](\t|$)! or next; m!\tC.pan:[^\s:]+(\t|$)! and next; print "$_\n"; }' | wc -l
# GrpSyn_000002   21      CApan:CaUC03G061460.1   CApan:CaUC03G061520.1   CApan:CaUC03G061550.1                   CApan:CaUC03G061560.1   CApan:Ciama_Chr03:24544226-24545011:+   CApan:Ciama_Chr03:30256016-30256508:+ 
# GrpSyn_000003   19      CApan:CaUC10G186420.1   CApan:CaUC10G186430.1   CApan:Ciama_Chr10:12115704-12116238:+   CLpan:Cla97C10G191540.1 CLpan:Cla97C10G191550.1 CLpan:Cla97C10G191553.1 CLpan:Cla97C10G191557.1

