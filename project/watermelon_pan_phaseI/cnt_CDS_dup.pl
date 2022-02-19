#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 in.bn6 > in.bn6.ident_filteredNum\n";

my @lvl = reverse(0 .. 100);

my %h;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  my $glob_matN = $ta[2]*$ta[3];
  for my $l1 (@lvl) {
    $glob_matN >= $ta[12] * $l1 or next;
    $h{$l1}{$ta[0]} ++;
  }
}
for my $l1 (@lvl) {
  $h{$l1} //= {};
  my $cnt = scalar(keys %{$h{$l1}});
  print "$l1\t$cnt\n";
}

# 
# [Sunhh@panda rmTE]$ head -3 1_maker_novCleanRmTEcomplete.c.fa.toRef.bn6 
# snap-NODE_38989__1_2429_ext-processed-gene-0.1-mRNA-1   Cla97C10G191640.1       93.83   81      3       2       439     517     295     375     5e-27   121     573     375     plus
# snap-NODE_38989__1_2429_ext-processed-gene-0.1-mRNA-1   Cla97C07G134700.1       93.83   81      3       2       439     517     295     375     5e-27   121     573     375     plus
# snap-NODE_38989__1_2429_ext-processed-gene-0.1-mRNA-1   Cla97C07G133870.1       94.52   73      2       2       439     509     295     367     3e-24   111     573     462     plus
#
