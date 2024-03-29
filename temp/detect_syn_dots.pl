#!/usr/bin/perl
# [4/12/2022] Make it simple and see if it works. I'll revise this if it is worthy.
use strict;
use warnings;
use Statistics::Regression;
use fileSunhh;

!@ARGV and die "perl $0 t1 out_prefix\n";

my $fn_xy = shift;
my $opref = shift;

# Load (x,y) points.
my @v_xy = &fileSunhh::load_tabFile( $fn_xy );
@v_xy = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @v_xy;
my @obs;
for (@v_xy) {
  push(@obs, [$_->[1], { 'const' => 1.0, 'someX' => $_->[0], 'ignored' => 'anything else' }  ]);
}
my %dist_p2p;
my %used;

# Set up cutoffs.
my $d3x = 300e3; # Maximum X distance of two neighboring points.
my $d3y = 300e3; # Maximum Y distance of two neighboring points.
my $d1 = sqrt($d3x**2 + $d3y**2); # Maximum point-to-point distance.
# $d1 = 500e3;
warn("d1=$d1\n");
my $d2 = 200e3; # Maximum point-to-line distance.
my $min_pts = 5;
my $tail_lmN = 5;

# Look for alignments from the first point.
my $blkN = 0;
open O1,'>',"$opref.sep" or die;
open O2,'>',"$opref.jn" or die;
for (my $i=0; $i<@v_xy; $i++) {
  defined $used{$i} and next;
  my @aln = ([$v_xy[$i][0], $v_xy[$i][1], $i]); # Currently, I allow one line for each point.
  &extAln(\@aln, \@v_xy, \%dist_p2p, \%used);
  scalar(@aln) == 0 and next;
  $blkN ++;
  my @segs;
  my ($min1, $min2, $max1, $max2);
  for my $a1 (@aln) {
    # print join("\t", $blkN, $a1->[0], $a1->[1], $aln[0][0], $aln[0][1], $aln[-1][0], $aln[-1][1])."\n";
    print O1 join("\t", $blkN, $a1->[0], $a1->[1])."\n";
    push(@segs, "$a1->[0],$a1->[1]");
    $min1 //= $a1->[0]; $min1 > $a1->[0] and $min1 = $a1->[0];
    $min2 //= $a1->[1]; $min2 > $a1->[1] and $min2 = $a1->[1];
    $max1 //= $a1->[0]; $max1 < $a1->[0] and $max1 = $a1->[0];
    $max2 //= $a1->[1]; $max2 < $a1->[1] and $max2 = $a1->[1];
  }
  print O2 join("\t", $blkN, scalar(@segs), $min1, $max1, $min2, $max2)."\n";
}
close O1;
close O2;


sub extAln {
  my ($aR, $vR, $dR, $usedR) = @_;
  # Look for the next available point which is nearest to the rightmost aligned point.
  my ($has_fit, $fA, $fB);
  $has_fit = 0;
  my $has_ext = 0;
  my @cand_pts;
  my $len_aR = scalar(@$aR);
  my $i = $aR->[-1][2];
  for (my $j=$i+1; $j<@$vR; $j++) {
    defined $usedR->{$j} and next;
    $vR->[$j][0] - $vR->[$i][0] > $d3x and last; # because this is sorted by X coords.
    abs($vR->[$j][1] - $vR->[$i][1]) > $d3y and next;
    unless (defined $dR->{$i}{$j}) {
      $dR->{$i}{$j} = &d_p2p( $vR->[$i], $vR->[$j] );
      $dR->{$j}{$i} = $dR->{$i}{$j};
    }
    $dR->{$i}{$j} > $d1 and next;
    if ($len_aR == 1) {
      push(@cand_pts, [$dR->{$i}{$j}, $j]);
      $has_ext = 1;
      next;
    }
    if ($has_fit == 0) {
      my $reg = Statistics::Regression->new( "pain", [ "const", "someX"] );
      for (my $ti = $#$aR; $ti > $#$aR-$tail_lmN and $ti >= 0; $ti--) {
        $reg->include( $obs[$aR->[$ti][2]][0], $obs[$aR->[$ti][2]][1] );
      }
      ($fA, $fB) = $reg->theta();
      $has_fit = 1;
    }
    my $t_d2 = &d_p2l($vR->[$j], $fA, $fB);
    $t_d2 > $d2 and next;
    # Choose the one with the lowest point-to-line distance!
    # push(@cand_pts, [$t_d2, $j]);
    push(@cand_pts, [$dR->{$i}{$j}, $j]);
    $has_ext = 1;
  }
  if ($has_ext == 1) {
    @cand_pts = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @cand_pts;
    my $j = $cand_pts[0][1];
    push(@$aR, [$vR->[$j][0], $vR->[$j][1], $j]);
    &extAln($aR, $vR, $dR, $usedR);
  }

  # Check if this extended alignment is good.
  if ( scalar(@$aR) >= $min_pts ) {
    for (@$aR) {
      $usedR->{$_->[2]} = 1;
    }
  } else {
    @$aR = ();
  }
  
  return;
}# extAln()

sub d_p2l {
  my ($p1, $a, $b) = @_;
  my $numerator = $p1->[1] - ($a + $b * $p1->[0]);
  my $denominator = sqrt(1 + $b**2);
  return( $numerator / $denominator );
}# d_p2l()

sub d_p2p {
  my ($p1, $p2) = @_;
  return(sqrt(($p1->[0]-$p2->[0])**2+($p1->[1]-$p2->[1])**2));
}# d_p2p()

# # Add data points
# $reg->include( 2.0, [ 1.0, 3.0, -1.0 ] );
# $reg->include( 1.0, [ 1.0, 5.0, 2.0 ] );
# $reg->include( 20.0, [ 1.0, 31.0, 0.0 ] );
# $reg->include( 15.0, [ 1.0, 11.0, 2.0 ] );
# 
# my %d;
# $d{const} = 1.0; $d{someX}= 5.0; $d{someY}= 2.0; $d{ignored}="anything else";
# $reg->include( 3.0, \%d );  # names are picked off the Regression specification
# 
# # Finally, print the result
# $reg->print();


