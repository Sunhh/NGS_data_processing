#!/usr/bin/perl
# 6/2/2023: If two SAM alignments neighbor each other on both query and reference genome, join them.
#   The input SAM file should have been sorted by position! No softclipping (S) is allowed!

use strict;
use warnings;
use fileSunhh;
use SeqAlnSunhh;

-t and !@ARGV and die "perl $0 in_seg.sam > joined.sam\n";

my $rev_flag = &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=1');
my $fwd_flag = &SeqAlnSunhh::mk_flag('keep'=>'2=0,4=0');
my %flag2str;
for (keys %$rev_flag) {$flag2str{$_}='r';}
for (keys %$fwd_flag) {$flag2str{$_}='f';}

my (%p);
while (<>) {
  chomp;
  if (m!^\@!) {
    print STDOUT "$_\n";
    next;
  }
  my @ta=split(/\t/, $_);
  $ta[2] eq '*' and next;
  defined $flag2str{$ta[1]} or die "bad 4: flag: $ta[1]\n";
  my $currStr   = $flag2str{$ta[1]};
  $ta[5] =~ m!S! and die "bad 3: No softclipping (S) is accepted!\n";
  my $cR = &SeqAlnSunhh::parseCigar($ta[5]);
  my $currRS = $ta[3];
  my $currRE = $ta[3] + $cR->{'SpanRefLen'} - 1;
  my ($currQS, $currQE);
  if ($ta[5] =~ m!^(\d+)H!) {
    $currQS = $1+1;
  } else {
    $currQS = 1;
  }
  if ($ta[5] =~ m!(\d+)H$!) {
    $currQE = $cR->{'RdLen'} - $1;
  } else {
    $currQE = $cR->{'RdLen'};
  }
  if ($currStr eq 'r') {
    my $s1 = $cR->{'RdLen'}-$currQS+1;
    my $e1 = $cR->{'RdLen'}-$currQE+1;
    ($currQE, $currQS) = ($s1, $e1); # Always QS <= QE;
  }
  if (defined $p{'RID'} and $p{'RID'} eq $ta[2] and $p{'QID'} eq $ta[0] and $p{'str'} eq 'f' and $currStr eq 'f' and $p{'RE'} == $currRS-1 and $p{'QE'} == $currQS-1 and $ta[9] ne '*' and $p{'line'}[9] ne '*') {
    $p{'QE'} = $currQE;
    $p{'RE'} = $currRE;
    $p{'line'}[5] =~ s!\d+H$!!;
    $ta[5] =~ s!^\d+H!!;
    $p{'line'}[5] =~ s!(\d+)([A-Z])$!! or die "bad 1: [$p{'line'}[5]]\n";
    my ($v1,$c1) = ($1,$2);
    if ($ta[5] =~ s!^(\d+)$c1!!) {
      $p{'line'}[5] .= ($v1+$1).$c1.$ta[5];
    } else {
      $p{'line'}[5] .= "$v1$c1$ta[5]";
    }
    $p{'line'}[4] < $ta[4] and $p{'line'}[4] = $ta[4];
    $p{'line'}[9] .= $ta[9];
    for my $tb (@ta[11 .. $#ta]) {
      if ($tb =~ m!^NM:i:(\d+)$!) {
        my $nm = $1;
        for my $tc (@{$p{'line'}}[11 .. $#{$p{'line'}}]) {
          if ($tc =~ m!^NM:i:(\d+)$!) {
            $tc = "NM:i:".($1+$nm);
            last;
          }
        }
        last;
      }
    }
  } elsif (defined $p{'RID'} and $p{'RID'} eq $ta[2] and $p{'QID'} eq $ta[0] and $p{'str'} eq 'r' and $currStr eq 'r' and $p{'RE'} == $currRS-1 and $p{'QS'} == $currQE+1 and $ta[9] ne '*' and $p{'line'}[9] ne '*') {
    $p{'RE'} = $currRE;
    $p{'QS'} = $currQS;
    $p{'line'}[5] =~ s!\d+H$!!;
    $ta[5] =~ s!^\d+H!!;
    $p{'line'}[5] =~ s!(\d+)([A-Z])$!! or die "bad 1a: [$p{'line'}[5]]\n";
    my ($v1,$c1) = ($1,$2);
    if ($ta[5] =~ s!^(\d+)$c1!!) {
      $p{'line'}[5] .= ($v1+$1).$c1.$ta[5];
    } else {
      $p{'line'}[5] .= "$v1$c1$ta[5]";
    }
    $p{'line'}[4] < $ta[4] and $p{'line'}[4] = $ta[4];
    $p{'line'}[9] .= $ta[9];
    for my $tb (@ta[11 .. $#ta]) {
      if ($tb =~ m!^NM:i:(\d+)$!) {
        my $nm = $1;
        for my $tc (@{$p{'line'}}[11 .. $#{$p{'line'}}]) {
          if ($tc =~ m!^NM:i:(\d+)$!) {
            $tc = "NM:i:".($1+$nm);
            last;
          }
        }
        last;
      }
    }
  } else {
    defined $p{'RID'} and print join("\t", @{$p{'line'}})."\n";
    %p=();
    $p{'RID'} = $ta[2];
    $p{'QID'} = $ta[0];
    $p{'str'} = $currStr;
    $p{'RS'}  = $currRS;
    $p{'RE'}  = $currRE;
    $p{'QS'}  = $currQS;
    $p{'QE'}  = $currQE;
    $p{'line'} = [@ta];
  }
}
if (defined $p{'RID'}) {
  print join("\t", @{$p{'line'}})."\n";
  %p=();
}

