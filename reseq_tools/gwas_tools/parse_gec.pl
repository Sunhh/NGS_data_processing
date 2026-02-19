#!/usr/bin/perl
use strict;
use warnings;

my $htxt = <<HH;
############################################################
perl $0 snp_gec.sum > snp_gec.sum.cuts
After CMD: java -jar gecV0.2/gec.jar --effect-number --plink-binary snp --genome --out snp_gec
############################################################
HH

-t and !@ARGV and die $htxt;

my ($cut1,$cut2) = (-1,-1);
while (<>) {
  chomp; $. == 1 and next; my @a=split(/\t/, $_);
  $cut1 = -1*log($a[4])/log(10);
  $cut2 = -1*log($a[3])/log(10);
  last;
}
print STDOUT join("\t", qw/cut1_signif cut2_suggest/)."\n";
print STDOUT join("\t", $cut1, $cut2)."\n";
# &fileSunhh::write2file("$dirN/snp.cuts", "$cut1\t$cut2\n", ">");

