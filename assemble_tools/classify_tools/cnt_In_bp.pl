#!/usr/bin/perl
use strict; 
use warnings; 


my %inCnt;
my %ord;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'qseqid' and next;
  $ta[6] =~ m!In:([\.\d]+)! or next;
  my $inDep = $1;
  $inDep > 0 or next;
  $inCnt{$ta[0]} += $ta[4];
  $ord{$ta[0]} //= $.;
}
for my $k1 (sort { $ord{$a} <=> $ord{$b} } keys %ord) {
  print join("\t", $)."\n";
}



# Sunhh@swift:/data/Sunhh/cmaxima/rmcont$ head -4 hf2.noRed.tochk.fa.toNt.bn6.jnInEx
#qseqid  qlen    qstart  qend    qspan   KingdomCounts   InExcludeCounts
#ptg000012l      137195  1       7424    7424    rDNA:5.28       Ex:5.28
#ptg000012l      137195  7425    7472    48      Un:0.00 In:0.00
#ptg000012l      137195  7473    18297   10825   rDNA:16.51      Ex:16.51

