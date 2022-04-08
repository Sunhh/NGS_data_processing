#!/usr/bin/perl
# [3/26/2022] There are genes/mRNAs with exons on different chromosomes. I want to choose only one with longest CDS as representative in cross-species analysis.
use strict;
use warnings;

-t and !@ARGV and die "perl $0 WCGv2.chr.gff3 > WCGv2.chr.trim2CDS.gff3\n";

my %rec;
my $mID='NA';
my %recG;
while (<>) {
  m!^\s*#|^\s*$! and next;
  chomp;
  my @ta=split(/\t/, $_);
  if      ($ta[2] =~ m!^gene$!i) {
    $ta[8] =~ m!^ID=([^\s;]+)! or die "$ta[8]\n";
    my $gID = $1;
    push(@{$recG{$gID}}, [@ta]);
  } elsif ($ta[2] =~ m!^mRNA$!i) {
    $ta[8] =~ m!^ID=([^\s;]+)! or die "$ta[8]\n";
    $mID = $1;
    $rec{$mID}{'rank'} //= $.;
    defined $rec{$mID}{'mLine'} and die "repeat mID $mID\n";
    push(@{$rec{$mID}{'mLine'}}, [@ta]);
    push(@{$rec{$mID}{'str'}}, $ta[6]);
    $rec{$mID}{'cLine'} //= [];
    $ta[8] =~ m!(?:^|;\s*)Parent=([^\s;]+)!i or die "mrna:$ta[8]\n";
    $rec{$mID}{'gID'} = $1;
  } elsif ($ta[2] =~ m!^CDS$!i) {
    $ta[8] =~ m!(?:^|;\s*)Parent=([^\s;]+)!i or die "cds:$ta[8]\n";
    my $c_mID = $1;
    $rec{$c_mID}{'cLine'} //= [];
    $rec{$c_mID}{'rank'} //= $.;
    push(@{$rec{$c_mID}{'cLine'}}, [@ta]);
    # $rec{$c_mID}{'cdsS'} //= $ta[3]; $rec{$c_mID}{'cdsS'} > $ta[3] and $rec{$c_mID}{'cdsS'} = $ta[3];
    # $rec{$c_mID}{'cdsE'} //= $ta[4]; $rec{$c_mID}{'cdsE'} < $ta[4] and $rec{$c_mID}{'cdsE'} = $ta[4];
  }
}
for my $mID (sort {$rec{$a}{'rank'} <=> $rec{$b}{'rank'}} keys %rec) {
  scalar(@{$rec{$mID}{'cLine'}}) > 0 or next;
  # Count how many chromosomes we need
  my %chrInfo; # {chrID} = [start, end, strand, rank]
  {
    my $cCnt = 0;
    for my $cL (@{$rec{$mID}{'cLine'}}) {
      $cCnt++;
      $chrInfo{$cL->[0]} //= [@{$cL}[3,4,6], $cCnt, $cL->[4]-$cL->[3]+1]; # [start, end, strand, rank]
      $chrInfo{$cL->[0]}[0] > $cL->[3] and $chrInfo{$cL->[0]}[0] = $cL->[3];
      $chrInfo{$cL->[0]}[1] < $cL->[4] and $chrInfo{$cL->[0]}[1] = $cL->[4];
      $chrInfo{$cL->[0]}[2] eq $cL->[6] or die "[Err] Bad ID: $mID\n";
      $chrInfo{$cL->[0]}[4] += ($cL->[4]-$cL->[3]+1);
    }
  }
  my ($best_chrID) = (sort {$chrInfo{$b}[4] <=> $chrInfo{$a}[4]} keys %chrInfo);
  # Retrieve general gene ID;
  defined $rec{$mID}{'gID'} or $rec{$mID}{'gID'} = "${mID}_hsG";
  my $gID = $rec{$mID}{'gID'};
  # Fill gene and mRNA lines for the best (chosen) chromosome.
  my (@gLine, @mLine);
  {
    for my $gL (@{$recG{$gID}}) {
      $gL->[0] eq $best_chrID or next;
      @gLine = @$gL;
      last;
    }
    scalar(@gLine) == 0 and do { @gLine = ($best_chrID, "fit", "gene", $chrInfo{$best_chrID}[0], $chrInfo{$best_chrID}[1], ".", $chrInfo{$best_chrID}[2], ".", "ID=$gID"); };
    for my $mL (@{$rec{$mID}{'mLine'}}) {
      $mL->[0] eq $best_chrID or next;
      @mLine = @$mL;
      last;
    }
    scalar(@mLine) == 0 and do { @mLine = ($best_chrID, "fit", "mRNA", $chrInfo{$best_chrID}[0], $chrInfo{$best_chrID}[1], ".", $chrInfo{$best_chrID}[2], ".", "ID=$mID;Parent=$gID"); };
    @gLine[3,4] = ($chrInfo{$best_chrID}[0], $chrInfo{$best_chrID}[1]);
    @mLine[3,4] = ($chrInfo{$best_chrID}[0], $chrInfo{$best_chrID}[1]);
  }
  print STDOUT join("\t", @gLine)."\n";
  print STDOUT join("\t", @mLine)."\n";
  for my $cL (@{$rec{$mID}{'cLine'}}) {
    $cL->[0] eq $best_chrID or next;
    print STDOUT join("\t", @$cL)."\n";
  }
}

