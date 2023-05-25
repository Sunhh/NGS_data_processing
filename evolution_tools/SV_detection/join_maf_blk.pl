#!/usr/bin/perl
### The alignments should be non-overlapping!!!
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 align2w38_anc_fix.good.maf > align2w38_anc_jn.maf\n";

my %aln;

while (<>) {
  chomp;
  if (m!^\s*#!) {
    print STDOUT "$_\n";
    next;
  }
  my $l2=<>; chomp($l2);
  my $l3=<>; chomp($l3);
  my $l4=<>; chomp($l4);
  m!^a\s! or die "bad2:$_\n";
  $l2 =~ m!^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)$! or die "bad3:$l2\n";
  my ($rID, $rS, $rLen, $rStr, $rSize, $rAln) = ($1,$2,$3,$4,$5,$6);
  $l3 =~ m!^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)! or die "bad4:$l3\n";
  my ($qID, $qS, $qLen, $qStr, $qSize, $qAln) = ($1,$2,$3,$4,$5,$6);
  push(@{$aln{$rID}{$qStr}}, [ [$rID, $rS, $rLen, $rSize, $rAln], [$qID, $qS, $qLen, $qSize, $qAln] ]);
}
# Sort %aln
for my $rID (sort keys %aln) {
  if (defined $aln{$rID}{'+'}) {
    @{$aln{$rID}{'+'}} = sort { $a->[1][0] cmp $b->[1][0] || $a->[0][1] <=> $b->[0][1] || $a->[1][1] <=> $b->[1][1] } @{$aln{$rID}{'+'}};
    my $tm = $aln{$rID}{'+'};
    my @prev = @{$tm->[0]};
    for (my $i=1; $i<@$tm; $i++) {
      my $tr = $tm->[$i];
      if ($prev[1][0] eq $tr->[1][0] and $prev[0][1]+$prev[0][2] == $tr->[0][1] and $prev[1][1]+$prev[1][2] == $tr->[1][1]) {
        # Same Q_ID, continuous R block, continuous Q block.(same strand)
        $prev[0][2] += $tr->[0][2];
        $prev[0][4] .= $tr->[0][4];
        $prev[1][2] += $tr->[1][2];
        $prev[1][4] .= $tr->[1][4];
      } else {
        # Output MAF alignment in prev.
        print STDOUT "a\tscore=0\n";
        print STDOUT join("\t", "s", $prev[0][0], $prev[0][1], $prev[0][2], "+", $prev[0][3], $prev[0][4])."\n";
        print STDOUT join("\t", "s", $prev[1][0], $prev[1][1], $prev[1][2], "+", $prev[1][3], $prev[1][4])."\n";
        print STDOUT "\n";
        # Renew @prev
        @prev = @$tr;
      }
    }
    # Output MAF alignment in @prev.
    print STDOUT "a\tscore=0\n";
    print STDOUT join("\t", "s", $prev[0][0], $prev[0][1], $prev[0][2], "+", $prev[0][3], $prev[0][4])."\n";
    print STDOUT join("\t", "s", $prev[1][0], $prev[1][1], $prev[1][2], "+", $prev[1][3], $prev[1][4])."\n";
    print STDOUT "\n";
    @prev=();
  }
  if (defined $aln{$rID}{'-'}) {
    @{$aln{$rID}{'-'}} = sort { $a->[1][0] cmp $b->[1][0] || $a->[0][1] <=> $b->[0][1] || $b->[1][1] <=> $a->[1][1] } @{$aln{$rID}{'-'}};
    my $tm = $aln{$rID}{'-'};
    my @prev=@{$tm->[0]};
    for (my $i=1; $i<@$tm; $i++) {
      my $tr = $tm->[$i];
      if ($prev[1][0] eq $tr->[1][0] and $prev[0][1]+$prev[0][2] == $tr->[0][1] and $prev[1][1]+$prev[1][2] == $tr->[1][1]) {
        $prev[0][2] += $tr->[0][2];
        $prev[0][4] .= $tr->[0][4];
        $prev[1][2] += $tr->[1][2];
        $prev[1][4] .= $tr->[1][4];
      } else {
        # Output MAF alignment in prev.
        print STDOUT "a\tscore=0\n";
        print STDOUT join("\t", "s", $prev[0][0], $prev[0][1], $prev[0][2], "+", $prev[0][3], $prev[0][4])."\n";
        print STDOUT join("\t", "s", $prev[1][0], $prev[1][1], $prev[1][2], "-", $prev[1][3], $prev[1][4])."\n";
        print STDOUT "\n";
        # Renew @prev
        @prev = @$tr;
      }
    }
    # Output MAF alignment in @prev.
    print STDOUT "a\tscore=0\n";
    print STDOUT join("\t", "s", $prev[0][0], $prev[0][1], $prev[0][2], "+", $prev[0][3], $prev[0][4])."\n";
    print STDOUT join("\t", "s", $prev[1][0], $prev[1][1], $prev[1][2], "-", $prev[1][3], $prev[1][4])."\n";
    print STDOUT "\n";
    @prev=();
  }
}

