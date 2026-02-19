#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 minuslogP_cut in_fastlmm.res peak_distance_1000000 > in_fastlmm.res.peak\n";

my $pCut=shift;
my $fn=shift;
my $dist=1e6;
scalar(@ARGV) > 0 and $dist = shift;

#my $rPcut = $pCut;
if (-f $pCut) {
  my @a1=&fileSunhh::load_tabFile($pCut);
  $a1[0][0] eq 'cut1_signif' or die "[Err] Bad input_pCut [$pCut]\n";
  $pCut = 10**(-$a1[1][0]);
} else {
  $pCut = 10**(-$pCut);
}
# warn "[Msg] cutoff=|$pCut|\n";
#die "[$pCut][$v][$rPcut]\n";
my $ifh=&openFH($fn);
my $ofh=&openFH("$fn.sig1", '>');
my @vars;
while (<$ifh>) {
  chomp;
  my @ta=split(/\t/,$_);
  # $ta[5] <= $pCut or next; # For emmax .p file.
  $ta[0] eq 'SNP' and next;
  $ta[4] <= $pCut or next; # For FaSTLMM 
  push(@vars, [@ta]);
  print {$ofh} "$_\n";
}
close($ofh);
close($ifh);

my @peaks;
for my $l (sort {$a->[4] <=> $b->[4]} @vars) {
  my $is=1;
  for my $pt (@peaks) {
    $pt->[1] eq $l->[1] or next; # chromosome ID; 1 in fastlmm, 0 in emmax.
    abs($pt->[3] - $l->[3])+1 > $dist and next; # position. Distance to peak. 3 in fastlmm, 1 in emmax.
    $is=0; last;
  }
  if ($is == 1) {
    push(@peaks, [@$l]);
    print STDOUT join("\t", @$l)."\n";
  }
}

