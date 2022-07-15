#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 min_var_len_20 C38.ndf_ref_struct.rmGap.effG > C38.ndf_ref_struct.rmGap.effG.slct1\n";

my $min_var_len = shift; # 20

while (<>) {
  m!^\s(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  m!\tVAR_annot|Name=(insertion|duplication|tandem_duplication|unaligned_end|unaligned_beginning|deletion|collapsed_repeat|collapsed_tandem_repeat);! or next;
  if (m!(?:overlap|subst|del|ins|blk)_len=(\d+)!){
    $1 >= $min_var_len or next;
  }
  print STDOUT "$_\n";
}

