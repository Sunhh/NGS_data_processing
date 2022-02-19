#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 file_list > summary.tbl\n";

my @lvl_IDPs;
for (my $i=1; $i<=5; $i++) {
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
    if (m!^(\S+)\t(\d+)\t(\d+)\t(\d+)$!) {
      my @ta=($1,$2,$3,$4);
      $fnCnt ++; $fnRank{$fn0} //= $fnCnt;
      &proc_line(\%h, $fn0, \@ta, \@lvl_IDPs);
    } else {
      # This file is a list file.
      m!^(\S+)! or die "$_\n";
      my $tafn = $1;
      my $fh1 = &openFH($tafn, '<');
      &tsmsg("[Msg] Processing file [$tafn]\n");
      $fnCnt ++; $fnRank{$tafn} ++;
      while (my $l1 = <$fh1>) {
        chomp($l1);
        $l1 =~ m!^(\S+)\t(\d+)\t(\d+)\t(\d+)$! or die "[Err] Bad line: $l1\n";
        my @tb = ($1,$2,$3,$4);
        &proc_line(\%h, $tafn, \@tb, \@lvl_IDPs);
      }
      close($fh1);
    }
  }
  close($fh0);
}
print STDOUT join("\t", 'filename', (map { "D${_}" } @lvl_IDPs))."\n";
for my $fn (sort {$fnRank{$a} <=> $fnRank{$b}} keys %h) {
  my @o = ($fn);
  my $ttlN;
  for my $lvlV (@lvl_IDPs) {
    $h{$fn}{$lvlV}{'P'} //= 0;
    push(@o, $h{$fn}{$lvlV}{'P'});
  }
  print STDOUT join("\t", @o)."\n";
}
sub proc_line {
  my ($hR, $fn, $lineR, $lvlR) = @_;
  for (my $i=0; $i<@$lvlR; $i++) {
    if ($lineR->[3] >= $lvlR->[$i]) {
      $hR->{$fn}{$lvlR->[$i]}{'P'} += $lineR->[2]-$lineR->[1];
    }
  }
}# proc_line()


# [Sunhh@panda 05_bwa_aln]$ head list.1 
# cnt_cov/ASM1001_CharlestonGray.CLpan.bam.gcovH
# cnt_cov/ASM1002_97103.CLpan.bam.gcovH
# cnt_cov/WM105_TS409MU.CLpan.bam.gcovH
# [Sunhh@panda 05_bwa_aln]$ head -4 cnt_cov/WM147_TuBaiPi.CLpan.bam.gcovH
# Cla97Chr01	11	43	1
# Cla97Chr01	43	58	2
# Cla97Chr01	58	66	3

