#!/usr/bin/perl
# [2/8/2022] Change the assembled contig locations from rmCont positions to dedup positions.
use strict; 
use warnings; 

!@ARGV and die "perl $0 Clean_UA.CD90.NR.fa.map2Src.tbl > map.dedup_rmcont_asm.tbl\n"; 

my $f1_tbl = shift; 

open F1,'<',"$f1_tbl" or die; 
while (<F1>) {
  m!^\s*$! and next;
  chomp; 
  my @ta=split(/\t/, $_); 
  if ($. == 1 and $ta[1] eq 'rmcont_ID') {
    print STDOUT join("\t", $_, qw/acc_ID asm_ID asm_start asm_end rmcont_startInAsm rmcont_endInAsm/)."\n"; 
    next; 
  }
  my ($accID, $asmID, $asmS, $asmE) = ('NA', 'NA', 'NA', 'NA'); 
  my ($ts, $te); 
  if      ( $ta[1] =~ m!^(\S+)_(NODE_\d+)F$!) {
    # Lei: Grif_15897_NODE_32941F ; WM1185_PI596668_NODE_10033F ;
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $ta[2], $ta[3]); 
  } elsif ( $ta[1] =~ m!^(\S+)_(NODE_\d+)__(\d+)_(\d+)F$! ) {
    # Lei/Honghe: 
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $3, $4); 
  } elsif ( $ta[1] =~ m!^((?:ASM|WM)\d+_[^\s_]+)_(\S+)__(\d+)_(\d+)F$! ) {
    # Lei/Honghe: (Not very possible) full sequence from a de novo assembly.
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $3, $4);
  } elsif ( $ta[1] =~ m!^((?:ASM|WM)\d+_[^\s_]+)_(\S+)F$! ) {
    # Lei: (Not very possible) full sequence from a de novo assembly.
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $ta[2], $ta[3]);
  } elsif ( $ta[1] =~ m!^(\S+)_(NODE_\d+)__(\d+)\-(\d+)[HMT]$! ) {
    # Lei: Grif_15897_NODE_1303__1-536H ; Grif_15897_NODE_18737__1703-2670T ; WM1185_PI596668_NODE_10044__637-1165M ; 
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $3, $4); 
  } elsif ( $ta[1] =~ m!^(\S+)_(NODE_\d+)__(\d+)_(\d+)__(\d+)\-(\d+)[HMT]$! ) {
    # Lei/Honghe: WM1185_PI596668_NODE_10044__1_2000__637-1165M; (may exist in the future.)
    ($accID, $asmID, $ts, $te, $asmS, $asmE) = ($1, $2, $3, $4, $5, $6); 
    $asmS = $ts + $asmS - 1;
    $asmE = $ts + $asmE - 1; 
  } elsif ( $ta[1] =~ m!^((?:ASM|WM)\d+_[^\s_]+)_(\S+)__(\d+)_(\d+)__(\d+)\-(\d+)[HMT]$! ) {
    # Lei/Honghe: ASM1001_CharlestonGray_CG_Chr00__10377199_10380681__1-1399H
    ($accID, $asmID, $ts, $te, $asmS, $asmE) = ($1, $2, $3, $4, $5, $6);
    $asmS = $ts + $asmS - 1;
    $asmE = $ts + $asmE - 1; 
  } elsif ( $ta[1] =~ m!^((?:ASM|WM)\d+_[^\s_]+)_(\S+)__(\d+)\-(\d+)[HMT]$! ) {
    # Lei: ASM1001_CharlestonGray_CG_Chr00__1-1399H (fake)
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $3, $4);
  } elsif ( $ta[1] =~ m!^(\S+)_(NODE_\d+)_(\d+)\-(\d+)$! ) {
    # Shan: WM1042_PI255136_NODE_10013_6821-7401;
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $3, $4);
  } elsif ( $ta[1] =~ m!^((?:ASM|WM)\d+_[^\s_]+)_(\S+)_(\d+)\-(\d+)$! ) {
    # Shan: WM1191_PI296341_chr0_1197_1-632
    ($accID, $asmID, $asmS, $asmE) = ($1, $2, $3, $4);
  } else {
    die "Failed to parse rmcont_ID [$ta[1]]\n"; 
  }
  my $asmS_dedup = $asmS+$ta[2]-1;
  my $asmE_dedup = $asmS+$ta[3]-1;
  print STDOUT join("\t", $_, $accID, $asmID, $asmS_dedup, $asmE_dedup, $asmS, $asmE)."\n"; 
}
close F1;


