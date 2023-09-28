#!/usr/bin/perl
# [3/28/2022] Look for orthologous gene pairs or gene locations.
#   To define a good gene mapping:
#   (1) At least 90% coverage transferred by Liftoff.
#   There are several special things to consider for gene-to-gene pairs:
#   (1) At least 100 bp or 50% overlapping of the both CDSs.
#     [3/29/2022]
#     In version 1, I used >= 100 bp OR 50% overlapping of the shorter one. 
#     But now I change to > 50% of both CDSs because I want to build a firm relationship at the very beginning, then I can further add those shorter genes back if they are left as singleton.
#   (2) Don't allow one gene to pair with two or more genes.
#     It is possible that one gene may cover two other genes, while one of the two overlaps smaller with another one. However, I would ignore this problem.
#   (3) For x-to-x pairs: (x-to-x is not 1-to-1)
#     Use the one with longest overlapping.
# [3/29/2022] Add local identity requirements.
# [9/28/2023] (1) Add options to assign min_cov, min_localID, min_ovlR, and allow_rep.
#             (2) With '-allow_rep', gene mapping can be multi-to-multi.
#             (3) When multiple R genes are mapped to the same Q/location, the added ones are not forcely converted to locations any more.
#                 In the past, I did this conversion to get one-to-one gene-to-gene relationships without synteny information.
#             (4) Allow use empty $fromXX and $toXX to ignore them.
#             (5) When -either_ovlR is given, instead of asking both overlapping ratio >= min_ovlR, any ratio >= min_ovlR approves a good gene pair.

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "min_cov:f",     # 0.9; Minimum coverage ratio of original gene body transferred by Liftoff.
  "min_localID:f", # 0.9; =identity/coverage Minimum local identification between original gene and transferred gene.
  "min_ovlR:f",    # 0.5; Minimum overlap ratio to say gene structures from two genomes are mapped.
  "either_ovlR!",
  "allow_rep!",    # If given, allow multi-to-multi gene/location pairs.
  "help!",
);

my $htxt = <<HHH;
################################################################
perl $0 from_tag to_tag output/mapCDS.CCpan.to.CLpan.tbl
 -min_cov      [0.9]
 -min_localID  [0.9]
 -allow_rep    [Boolean]
 -min_ovlR     [0.5]
 -either_ovlR  [Boolean]
################################################################
HHH

scalar(@ARGV) == 3 or die "$htxt";
$opts{'help'} and die "$htxt";

my $fromXX = shift;
my $toXX   = shift;

$fromXX ne '' and $fromXX !~ m!:$! and $fromXX .= ":";
$toXX   ne '' and $toXX   !~ m!:$! and $toXX   .= ":";

$opts{'min_cov'}     //= 0.9; my $min_cov = $opts{'min_cov'};
$opts{'min_localID'} //= 0.9; my $min_localID = $opts{'min_localID'};
$opts{'min_ovlR'}    //= 0.5; my $min_ovlR   = $opts{'min_ovlR'};

my @o1; # [ [ 0:ori_R_ID, 1:Q_ID/loc, 2:R2Q_ovlLen, 3:trans_R_ID, 4:trans_R_chrID, 5:trans_R_chrS, 6:trans_R_chrE, 7:trans_R_chrStr, 8:identity, 9:coverage, 10:trans_R_size, 11:trans_R_exNum, 12:ori_R_size, 13:ori_R_exNum, 14:Q_size, 15:Q_exNum, 16:Q_ID_type, 17:trans_R_CDSblks], ... ]
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'trans_R_ID' and next;
  # Eligible transferring.
  # (1) Require coverage > 0.5 and local identity > 0.9;
  $ta[1] >= $min_cov or next;
  $ta[2]/$ta[1] >= $min_localID or next;
  # Eligible gene pair.
  if      ($ta[9] eq '.') {
    # Not overlapping any genes.
    push(@o1, ["${fromXX}$ta[6]", "${toXX}$ta[12]:$ta[13]-$ta[14]:$ta[15]", $ta[3], $ta[0], @ta[12,13,14,15], @ta[2,1], @ta[4,5,7,8,10,11], $ta[9], $ta[24]]);
  } elsif ( $ta[3] >= $min_ovlR * $ta[7] and $ta[3] >= $min_ovlR * $ta[10] ) {
    # This is a good gene pair. 
    push(@o1, ["${fromXX}$ta[6]", "${toXX}$ta[9]",                          $ta[3], $ta[0], @ta[12,13,14,15], @ta[2,1], @ta[4,5,7,8,10,11], $ta[9], $ta[24]]);
  } elsif ( $opts{'either_ovlR'} and ($ta[3] >= $min_ovlR * $ta[7] or $ta[3] >= $min_ovlR * $ta[10]) ) {
    push(@o1, ["${fromXX}$ta[6]", "${toXX}$ta[9]",                          $ta[3], $ta[0], @ta[12,13,14,15], @ta[2,1], @ta[4,5,7,8,10,11], $ta[9], $ta[24]]);
  } else {
    # The overlapped region is smaller than the threshold.
    push(@o1, ["${fromXX}$ta[6]", "${toXX}$ta[12]:$ta[13]-$ta[14]:$ta[15]", $ta[3], $ta[0], @ta[12,13,14,15], @ta[2,1], @ta[4,5,7,8,10,11], ".",    $ta[24]]);
  }
}
@o1 = sort {
  my $res = $b->[8] <=> $a->[8] || $b->[2] <=> $a->[2]; # High to low: global identity, R2Q_ovlLen;
  # if ($res == 0 and $a->[10] ne ".") {
  if ($res == 0) {
    # I think '$a->[10]' (trans_R_size) should never be '.'.
    # Further sorted by [trans_R_size-ori_R_size] and [trans_R_exNum-ori_R_exNum];
    $res =  abs($a->[10]-$a->[12]) <=> abs($b->[10]-$b->[12]) || abs($a->[11]-$a->[13]) <=> abs($b->[11]-$b->[13]);
  }
  if ($res == 0) {
    # Further sorted by trans_R_chrID, trans_R_chrS, and trans_R_chrE
    $res = $a->[4] cmp $b->[4] || $a->[5] <=> $b->[5] || $a->[6] <=> $b->[6];
  }
  $res;
  } @o1;
# my @o1; # [ [ 0:ori_R_ID, 1:Q_ID/loc, 2:R2Q_ovlLen, 3:trans_R_ID, 4:trans_R_chrID, 5:trans_R_chrS, 6:trans_R_chrE, 7:trans_R_chrStr, 8:identity, 9:coverage, 10:trans_R_size, 11:trans_R_exNum, 12:ori_R_size, 13:ori_R_exNum, 14:Q_size, 15:Q_exNum, 16:Q_ID_type, 17:trans_R_CDSblks], ... ]
my @o2;
if ($opts{'allow_rep'}) {
  @o2 = @o1;
} else {
  my (%hasR, %hasQ);
  for my $to1 (@o1) {
    defined $hasR{$to1->[0]} and next;
    defined $hasQ{$to1->[1]} and next;
    $hasR{$to1->[0]} = 1;
    $hasQ{$to1->[1]} = 1;
    push(@o2, [@$to1]);
  }
  # Add lost in which there are multiple R genes mapping exactly the same Q gene/location.
  # For these R genes, I just use their mapping locations.
  for my $to1 (@o1) {
    defined $hasR{$to1->[0]} and next;
    $hasR{$to1->[0]} = 1;
    push(@o2, [@$to1]);
    # $o2[-1][1] = "${toXX}$to1->[4]:$to1->[5]-$to1->[6]:$to1->[7]";
  }
}
# Output mapping;
print STDOUT join("\t", qw/ori_R_ID Q_ID  R2Q_ovlLen  trans_R_ID  trans_R_chrID  trans_R_chrS  trans_R_chrE  trans_R_chrStr  identity  coverage  trans_R_size  trans_R_exNum  ori_R_size  ori_R_exNum  Q_size  Q_exNum Q_ID_type trans_R_CDSblks/)."\n";
for my $to2 (@o2) {
  print STDOUT join("\t", @$to2)."\n";
}

