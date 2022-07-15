#!/usr/bin/perl
# [5/6/2022] Use bedtools to retrieve overlapping information.
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;

my $htxt = <<HH;
##########################################################################################
# Calculate overlapping regions and what bed_samples these regions are carried.
#   bedtools program is required!!!
#   LogInforSunhh and fileSunhh modules are required.

perl $0 in_1.bed in_2.bed in_3.bed ... > ovl.cnt.bed

# Output: chrID \\t start (0-based) \\t end \\t # of samples \\t sample1;sample2;... / none
### "# of samples" shows how many samples carry this region.
### The first several lines initialized with '#' tell the sample labels and their relative bed file names.

# Workflow to prepare input files.
### Input: in_1.bam; in_2.bam; in_3.bam;
### commands:
###   bedtools genomecov -ibam in_1.bam -bg -max 5 > in_1.bam.gcov
###   bedtools genomecov -ibam in_2.bam -bg -max 5 > in_2.bam.gcov
###   bedtools genomecov -ibam in_3.bam -bg -max 5 > in_3.bam.gcov
###   perl cnt_ovl_fromBeds.pl in_1.bam.gcov  in_2.bam.gcov  in_3.bam.gcov  > overlap.cnt.bed

##########################################################################################
HH

!@ARGV and &LogInforSunhh::usage($htxt);

my $wdir = &fileSunhh::new_tmp_dir('create' => 1);

my (%boundaries, @lines, %rank_cid, @all_cid);
for (my $i=0; $i<@ARGV; $i++) {
  open F,'<',"$ARGV[$i]" or die;
  while (<F>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[1]++;
    push(@{$lines[$i]{$ta[0]}}, [@ta[1,2]]);
    defined $rank_cid{$ta[0]} or push(@all_cid, $ta[0]);
    $rank_cid{$ta[0]} //= $.;
    unless ( defined $boundaries{'has'}{$ta[0]}{$ta[1]} ) {
      push(@{$boundaries{'arr'}{$ta[0]}}, $ta[1]);
      $boundaries{'has'}{$ta[0]}{$ta[1]} = 1;
    }
    unless ( defined $boundaries{'has'}{$ta[0]}{$ta[2]} ) {
      push(@{$boundaries{'arr'}{$ta[0]}}, $ta[2]);
      $boundaries{'has'}{$ta[0]}{$ta[2]} = 1;
    }
  }
  close F;
}
my (%blks);
open O1,'>',"$wdir/basic.bed" or die;
for my $cid (@all_cid) {
  @{$boundaries{'arr'}{$cid}} = sort { $a<=>$b } @{$boundaries{'arr'}{$cid}};
  push(@{$blks{$cid}}, [$boundaries{'arr'}{$cid}[0], $boundaries{'arr'}{$cid}[0]]);
  for (my $i=1; $i<@{$boundaries{'arr'}{$cid}}; $i++) {
    if ( $boundaries{'arr'}{$cid}[$i-1] < $boundaries{'arr'}{$cid}[$i]-1 ) {
      push(@{$blks{$cid}}, [$boundaries{'arr'}{$cid}[$i-1]+1, $boundaries{'arr'}{$cid}[$i]-1]);
    }
  }
  scalar(@{$boundaries{'arr'}{$cid}}) > 1 and push(@{$blks{$cid}}, [$boundaries{'arr'}{$cid}[-1], $boundaries{'arr'}{$cid}[-1]]);
  for my $t1 (@{$blks{$cid}}) {
    my $ts = $t1->[0]-1;
    my $te = $t1->[1];
    print O1 join("\t", $cid, $ts, $te)."\n";
  }
}
close O1;

for (my $i=0; $i<@ARGV; $i++) {
  &runCmd("bedtools intersect -a $wdir/basic.bed -b $ARGV[$i] > $wdir/o_$i.bed");
  my %h1;
  open F1,'<',"$wdir/o_$i.bed" or die;
  while (<F1>) {
    chomp;
    $h1{$_} = 1;
  }
  close F1;
  for my $cid (@all_cid) {
    for my $t1 (@{$blks{$cid}}) {
      my $tkey = join("\t", $cid, $t1->[0]-1, $t1->[1]);
      if (defined $h1{$tkey}) {
        push(@{$t1->[2]}, $i);
      } else {
        $t1->[2] //= [];
      }
    }
  }
}
for (my $i=0; $i<@ARGV; $i++) {
  print STDOUT "# $i: $ARGV[$i]\n";
}
for my $cid (@all_cid) {
  for my $t1 (@{$blks{$cid}}) {
    my $txt = join(";", @{$t1->[2]}); $txt eq '' and $txt = "none";
    print STDOUT join("\t", $cid, $t1->[0]-1, $t1->[1], scalar(@{$t1->[2]}), $txt)."\n";
  }
}
&fileSunhh::_rmtree($wdir);

