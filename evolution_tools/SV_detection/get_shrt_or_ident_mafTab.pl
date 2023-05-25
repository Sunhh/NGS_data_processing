#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 38000 pident_99 align2w38_anc_fix.maf.blasttab > align2w38_anc_fix.maf.blasttab.shrt\n";

my $maxLen   = shift;
my $minIdent = shift;

while (<>) {
  chomp;
  my @ta=split(/\t/,$_);
  if (abs($ta[7]-$ta[6])+1 < $maxLen and abs($ta[9]-$ta[8])+1 < $maxLen) {
    print "$_\n";
  } elsif ($ta[2] >= $minIdent) {
    print "$_\n";
  }
}

# [Sunhh@panda with_allNewR2_ori]$ deal_table.pl align2w38_anc_fix.maf.blasttab -col_head 
# 0       22CEXU11_Chr06
# 1       22CEXU43_Chr06
# 2       100.00
# 3       2385
# 4       0
# 5       0
# 6       20755838
# 7       20758222
# 8       19857222
# 9       19859606

