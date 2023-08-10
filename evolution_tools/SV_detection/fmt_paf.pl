#!/usr/bin/perl
use strict;
use warnings;

-t and !@ARGV and die "perl $0 minimap2.paf > minimap2.paf.tbl\n";

print STDOUT join("\t", qw/QID QLen QS QE Str SID SLen SS SE Match BlkLen MapQ QIdent SIdent Ident/)."\n";
while (<>) {
  chomp;
  my @ta=split(/\t/,$_);
  # $ta[11] >= 20 or next;
  my $ident1 = sprintf("%0.1f", $ta[9]/($ta[3]-$ta[2])*100);
  my $ident2 = sprintf("%0.1f", $ta[9]/($ta[8]-$ta[7])*100);
  my $ident3 = sprintf("%0.1f", $ta[9]/$ta[10]*100);
  print join("\t", @ta[0..11], $ident1, $ident2, $ident3)."\n";
}

### head -2 PAF/refU43_22CEXU11.paf | deal_table.pl -transpose | deal_table.pl -label_mark 0..99 | less -S
# 0       22CEXU11_Chr01  22CEXU11_Chr01
# 1       37955663        37955663
# 2       36419029        37276016
# 3       37054497        37835476
# 4       +               +
# 5       22CEXU43_Chr01  22CEXU43_Chr01
# 6       37591536        37591536
# 7       36074564        36928847
# 8       36707037        37490474
# 9       607399          528954
# 10      650544          584416
# 11      60              60
# 12      tp:A:P          tp:A:P
# 13      mm:i:9998       mm:i:7717
# 14      gn:i:33147      gn:i:47745
# 15      go:i:2313       go:i:1723
# 16      cg:Z:77=1X83=1X61=1X110=1X11=1X11=6I1=1X11=1X7=2X22=
# 17      cs:Z::77*GA:83*AG:61*CT:110*CT:11*AG:11+CTCTCT:1*CA:

