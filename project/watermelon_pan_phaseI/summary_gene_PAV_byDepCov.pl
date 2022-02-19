#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 file_list > summary.tbl\n";

my @lvl_IDPs;
for (my $i=10; $i<=100; $i+=5) {
  push(@lvl_IDPs, $i);
}

my %h;
my %fnRank;
my $fnCnt = 0;
for my $fn0 (@ARGV) {
  my $fh0 = &openFH($fn0, '<');
  &tsmsg("[Msg] Processing file [$fn0]\n");
  while (<$fh0>) {
    m!^\s*(#|$)! and next;
    chomp;
    my @ta=split(/\t/, $_);
    if (scalar(@ta) >= 5) {
      # This file is a ".bed_gcov.intersect.pav_rd1_cov90" file.
      $fnCnt ++; $fnRank{$fn0} //= $fnCnt;
      &proc_line(\%h, $fn0, \@ta, \@lvl_IDPs);
    } else {
      # This file is a list file.
      my $fh1 = &openFH($ta[0], '<');
      &tsmsg("[Msg] Processing file [$ta[0]]\n");
      $fnCnt ++; $fnRank{$ta[0]} ++;
      while (my $l1 = <$fh1>) {
        chomp($l1);
        my @tb = split(/\t/, $l1);
        &proc_line(\%h, $ta[0], \@tb, \@lvl_IDPs);
      }
      close($fh1);
    }
  }
  close($fh0);
}
print STDOUT join("\t", qw/filename ttl_geneN/, (map { "Pr.ID_$_" } @lvl_IDPs), (map { "Pr.ID_${_}_Perc" } @lvl_IDPs))."\n";
for my $fn (sort {$fnRank{$a} <=> $fnRank{$b}} keys %h) {
  my @o = ($fn, 0);
  my $ttlN;
  for my $lvlV (@lvl_IDPs) {
    $h{$fn}{$lvlV}{'P'} //= 0;
    $h{$fn}{$lvlV}{'A'} //= 0;
    $ttlN //= $h{$fn}{$lvlV}{'P'} + $h{$fn}{$lvlV}{'A'};
    push(@o, $h{$fn}{$lvlV}{'P'});
  }
  $o[1] = $ttlN;
  print STDOUT join("\t", @o, (map { sprintf("%.2f", $_/$o[1]*100) } @o[2..$#o]))."\n";
}
sub proc_line {
  my ($hR, $fn, $lineR, $lvlR) = @_;
  $lineR->[0] eq 'geneID' and return;
  for (my $i=0; $i<@$lvlR; $i++) {
    if ($lineR->[3] >= $lvlR->[$i]/100 * $lineR->[2]) {
      $hR->{$fn}{$lvlR->[$i]}{'P'} ++;
    } else {
      $hR->{$fn}{$lvlR->[$i]}{'A'} ++;
    }
  }
}# proc_line()


# [Sunhh@panda 05_bwa_aln]$ head list.1 
# cnt_cov/ASM1001_CharlestonGray.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/ASM1002_97103.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM105_TS409MU.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM1176_PI179881.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM135_PI381740.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM147_TuBaiPi.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM148_KaLaQiPa.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM149_DaHongZi.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM150_ALaKeZiWai.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM151_XiaoZiGua.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# [Sunhh@panda 05_bwa_aln]$ head -4 cnt_cov/WM147_TuBaiPi.CLpan.bam.gcovH.intersect.pav_rd1_cov90 
# geneID  PAV     CDS_len CDS_cov CDS_covP
# Cla97C01G000010.1       P       762     762     100.00
# Cla97C01G000020.1       P       318     318     100.00
# Cla97C01G000030.1       P       1521    1521    100.00
