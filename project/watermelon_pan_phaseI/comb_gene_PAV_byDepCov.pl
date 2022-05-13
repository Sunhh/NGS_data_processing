#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 CDScoveragePercent(0-100] file_list > summary.tbl\n";

my $lvl_IDP = shift; 


my %h;
my %fnRank;
my $fnCnt = 0;
my $geneCnt = 0; 
my %basicGene;
for my $fn0 (@ARGV) {
  my $fh0 = &openFH($fn0, '<');
  &tsmsg("[Msg] Processing file [$fn0]\n");
  while (<$fh0>) {
    m!^\s*(#|$)! and next;
    chomp;
    my @ta=split(/\t/, $_);
    if (scalar(@ta) >= 5) {
      # This file is a ".bed_gcov.intersect.pav_rd1_cov90" file.
      $ta[0] eq 'geneID' and next;
      $fnCnt ++; $fnRank{$fn0} //= $fnCnt;
      $geneCnt ++;
      $basicGene{$ta[0]} //= [$ta[2], $geneCnt];
      &proc_line_1(\%h, $fn0, \@ta, $lvl_IDP);
    } else {
      # This file is a list file.
      my $fh1 = &openFH($ta[0], '<');
      &tsmsg("[Msg] Processing file [$ta[0]]\n");
      $fnCnt ++; $fnRank{$ta[0]} ++;
      while (my $l1 = <$fh1>) {
        chomp($l1);
        my @tb = split(/\t/, $l1);
        $tb[0] eq 'geneID' and next;
        $geneCnt ++;
        $basicGene{$tb[0]} //= [$tb[2], $geneCnt];
        &proc_line_1(\%h, $ta[0], \@tb, $lvl_IDP);
      }
      close($fh1);
    }
  }
  close($fh0);
}

my @FNs   = sort { $fnRank{$a} <=> $fnRank{$b} } keys %fnRank;
# genome_coverage/ASM1002_97103.CLpan.bam.gcov.intersect.pav_rd1_cov50
my @oAccs;
for my $fn1 (@FNs) {
  my $a1 = $fn1;
  $a1 =~ s!^genome_coverage/!!;
  $a1 =~ s!\.\w\wpan.bam.gcov.intersect.pav_rd\d+_cov\d+$!!;
  push(@oAccs, $a1);
}
my @genes = sort { $basicGene{$a}[1] <=> $basicGene{$b}[1] } keys %basicGene;
print STDOUT join("\t", "geneID", "CDS_len", "Pr_CDSrat", @oAccs, @oAccs)."\n";
for my $geneID (@genes) {
  my @o = ($geneID, $basicGene{$geneID}[0], $lvl_IDP);
  for my $fn1 (@FNs) {
    push(@o, $h{$geneID}{$fn1}[0]);
  }
  for my $fn1 (@FNs) {
    push(@o, $h{$geneID}{$fn1}[1]);
  }
  print STDOUT join("\t", @o)."\n";
}

sub proc_line_1 {
  my ($hR, $fn, $lineR, $lvl) = @_;
  $lineR->[0] eq 'geneID' and return;
  if ($lineR->[3] >= $lvl/100 * $lineR->[2]) {
    $hR->{$lineR->[0]}{$fn} = ['P', $lineR->[3]];
  } else {
    $hR->{$lineR->[0]}{$fn} = ['A', $lineR->[3]];
  }
  return;
}# proc_line_1()


# [Sunhh@panda 05_bwa_aln]$ head list.1 
# cnt_cov/ASM1001_CharlestonGray.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/ASM1002_97103.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# cnt_cov/WM151_XiaoZiGua.CLpan.bam.gcovH.intersect.pav_rd1_cov90
# [Sunhh@panda 05_bwa_aln]$ head -4 cnt_cov/WM147_TuBaiPi.CLpan.bam.gcovH.intersect.pav_rd1_cov90 
# geneID  PAV     CDS_len CDS_cov CDS_covP
# Cla97C01G000010.1       P       762     762     100.00
# Cla97C01G000020.1       P       318     318     100.00
# Cla97C01G000030.1       P       1521    1521    100.00
