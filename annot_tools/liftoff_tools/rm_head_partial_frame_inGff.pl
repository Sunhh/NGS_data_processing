#!/usr/bin/perl
# [01/21/2022] To make it easy, 
#   I don't accept multiple mRNAs in a same gene! 
#   And the element types must be pre-defined.
#   I don't accept non-linear gene models, which means the real exons must be mono ascending or mono decending, so the order of exons/CDS placement in the file does not matter with the structure.
# [5/25/2022] Delete the first CDS block if the initial base of the first good codon is in the second CDS block. Don't allow 5'-UTR if the first CDS block doesn't start with an intact codon.
use strict;
use warnings;
use mathSunhh;
my $ms_obj = mathSunhh->new();

!@ARGV and die "perl $0 in.gff3 > in_fix.gff3\n";

my $f2 = shift;

my %good_ele = qw(
gene  1
mrna  2
exon  3
five_prime_utr  4
cds             5
three_prime_utr 6
);

my ($f2_gHR, $f2_mHR) = &load_gff3($f2);
my %g2 = %$f2_gHR;
my %m2 = %$f2_mHR;

my %has_mrna;
# Step 1: Revise elements following gene lines.
for my $gID (sort { $g2{$a}{'lineN'}  <=> $g2{$b}{'lineN'} } keys %g2) {
  unless ( defined $g2{$gID}{'mrnaIDs'} ) {
    print STDOUT join("\t", @{$g2{$gID}{'geneLine'}})."\n";
    next;
  }
  my @o1;
  push(@o1, [@{$g2{$gID}{'geneLine'}}]);
  for my $mID (sort { $m2{$a}{'lineN'} <=> $m2{$b}{'lineN'} } keys %{$g2{$gID}{'mrnaIDs'}}) {
    # print STDOUT join("\t", @{$m2{$mID}{'mrnaLine'}})."\n";
    $has_mrna{$mID} = 1;
    $m2{$mID}{'cdsLines'} //= [];
    $m2{$mID}{'exonLines'} //= [];
    if (scalar(@{$m2{$mID}{'cdsLines'}}) > 0 and scalar(@{$m2{$mID}{'exonLines'}}) == 0) {
      for my $lR2 (@{$m2{$mID}{'cdsLines'}}) {
        push(@{$m2{$mID}{'exonLines'}}, [@$lR2]);
        $m2{$mID}{'exonLines'}[-1][2] = 'exon';
      }
    }
    $m2{$mID}{'5utrLines'} //= [];
    $m2{$mID}{'3utrLines'} //= [];
    CHK_G_CDSLINES:
    if (scalar(@{$m2{$mID}{'cdsLines'}}) > 0) {
      @{$m2{$mID}{'cdsLines'}} = sort { $a->[3] <=> $b->[3] } @{$m2{$mID}{'cdsLines'}};
      @{$m2{$mID}{'exonLines'}} = sort { $a->[3] <=> $b->[3] } @{$m2{$mID}{'exonLines'}};
      if      ($m2{$mID}{'cdsLines'}[-1][6] eq '-' and $m2{$mID}{'cdsLines'}[-1][7] =~ m!^[^.0]$!) {
        # It is not reasonable to have the 5'-UTR region and non-zero frame CDS initials, so it should be good to directly remove the bases in the exons.
        scalar(@{$m2{$mID}{'5utrLines'}}) > 0 and die "[Err] FFFF 8: [$mID] Should not have 5'-utr when the first codon is not intact.\n";
        $m2{$mID}{'exonLines'}[-1][4] == $m2{$mID}{'cdsLines'}[-1][4] or die "[Err] FFFF 7: [$mID] additional 5'-exons.\n";
        $m2{$mID}{'cdsLines'}[-1][4] -= $m2{$mID}{'cdsLines'}[-1][7];
        $m2{$mID}{'cdsLines'}[-1][7] = 0;
        $m2{$mID}{'exonLines'}[-1][4] = $m2{$mID}{'cdsLines'}[-1][4];
        $o1[0][4]                = $m2{$mID}{'cdsLines'}[-1][4];
        $m2{$mID}{'mrnaLine'}[4] = $m2{$mID}{'cdsLines'}[-1][4];
        if ($m2{$mID}{'cdsLines'}[-1][4] < $m2{$mID}{'cdsLines'}[-1][3]) {
          pop(@{$m2{$mID}{'cdsLines'}});
          pop(@{$m2{$mID}{'exonLines'}});
          $o1[0][4]                = $m2{$mID}{'cdsLines'}[-1][4];
          $m2{$mID}{'mrnaLine'}[4] = $m2{$mID}{'cdsLines'}[-1][4];
          goto CHK_G_CDSLINES;
        }
      } elsif ($m2{$mID}{'cdsLines'}[0][6] eq '+' and $m2{$mID}{'cdsLines'}[0][7]  =~ m!^[^.0]$!) {
        scalar(@{$m2{$mID}{'5utrLines'}}) > 0 and die "[Err] FFFF 9: [$mID] Should not have 5'-utr when the first codon is not intact.\n";
        $m2{$mID}{'exonLines'}[0][3] == $m2{$mID}{'cdsLines'}[0][3] or die "[Err] FFFF 7p: [$mID] additional 5'-exons.\n";
        $m2{$mID}{'cdsLines'}[0][3] += $m2{$mID}{'cdsLines'}[0][7];
        $m2{$mID}{'cdsLines'}[0][7] = 0;
        $m2{$mID}{'exonLines'}[0][3] = $m2{$mID}{'cdsLines'}[0][3];
        $o1[0][3]                = $m2{$mID}{'cdsLines'}[0][3];
        $m2{$mID}{'mrnaLine'}[3] = $m2{$mID}{'cdsLines'}[0][3];
        if ($m2{$mID}{'cdsLines'}[0][3] > $m2{$mID}{'cdsLines'}[0][4]) {
          # This CDS block should be removed entirely.
          shift(@{$m2{$mID}{'cdsLines'}});
          shift(@{$m2{$mID}{'exonLines'}});
          $o1[0][3]                = $m2{$mID}{'cdsLines'}[0][3];
          $m2{$mID}{'mrnaLine'}[3] = $m2{$mID}{'cdsLines'}[0][3];
          goto CHK_G_CDSLINES;
        }
      }
    }
    my @o2;
    if ($m2{$mID}{'mrnaLine'}[6] eq '-') {
      @o2 = sort { $b->[4] <=> $a->[4] || $good_ele{lc($a->[2])} <=> $good_ele{lc($b->[2])} } (
        @{$m2{$mID}{'exonLines'}},@{$m2{$mID}{'5utrLines'}},@{$m2{$mID}{'cdsLines'}},@{$m2{$mID}{'3utrLines'}}
      );
    } else {
      @o2 = sort { $a->[3] <=> $b->[3] || $good_ele{lc($a->[2])} <=> $good_ele{lc($b->[2])} } (
        @{$m2{$mID}{'exonLines'}},@{$m2{$mID}{'5utrLines'}},@{$m2{$mID}{'cdsLines'}},@{$m2{$mID}{'3utrLines'}}
      );
    }
    push(@o1, [@{$m2{$mID}{'mrnaLine'}}], @o2);
  }
  for my $lR2 (@o1) {
    print STDOUT join("\t", @$lR2)."\n";
  }
}

# Step 2: Revise elements following only mRNA lines without a gene line (parent).
for my $mID (sort { $m2{$a}{'lineN'} <=> $m2{$b}{'lineN'} } keys %m2) {
  defined $has_mrna{$mID} and next;
  $m2{$mID}{'cdsLines'} //= [];
  $m2{$mID}{'exonLines'} //= [];
  if (scalar(@{$m2{$mID}{'cdsLines'}}) > 0 and scalar(@{$m2{$mID}{'exonLines'}}) == 0) {
    for my $lR2 (@{$m2{$mID}{'cdsLines'}}) {
      push(@{$m2{$mID}{'exonLines'}}, [@$lR2]);
      $m2{$mID}{'exonLines'}[-1][2] = 'exon';
    }
  }
  $m2{$mID}{'5utrLines'} //= [];
  $m2{$mID}{'3utrLines'} //= [];
  CHK_M_CDSLINES:
  if (scalar(@{$m2{$mID}{'cdsLines'}}) > 0) {
    @{$m2{$mID}{'cdsLines'}} = sort { $a->[3] <=> $b->[3] } @{$m2{$mID}{'cdsLines'}};
    @{$m2{$mID}{'exonLines'}} = sort { $a->[3] <=> $b->[3] } @{$m2{$mID}{'exonLines'}};
    if      ($m2{$mID}{'cdsLines'}[-1][6] eq '-' and $m2{$mID}{'cdsLines'}[-1][7] =~ m!^[^.0]$!) {
      scalar(@{$m2{$mID}{'5utrLines'}}) > 0 and die "[Err] FFFF 8Mm: [$mID] Should not have 5'-utr when the first codon is not intact.\n";
      $m2{$mID}{'exonLines'}[-1][4] == $m2{$mID}{'cdsLines'}[-1][4] or die "[Err] FFFF 7Mm: [$mID] additional 5'-exons.\n";
      $m2{$mID}{'cdsLines'}[-1][4] -= $m2{$mID}{'cdsLines'}[-1][7];
      $m2{$mID}{'cdsLines'}[-1][7] = 0;
      $m2{$mID}{'exonLines'}[-1][4] = $m2{$mID}{'cdsLines'}[-1][4];
      $m2{$mID}{'mrnaLine'}[4]      = $m2{$mID}{'cdsLines'}[-1][4];
      if ($m2{$mID}{'cdsLines'}[-1][4] < $m2{$mID}{'cdsLines'}[-1][3]) {
        pop(@{$m2{$mID}{'cdsLines'}});
        pop(@{$m2{$mID}{'exonLines'}});
        $m2{$mID}{'mrnaLine'}[4]      = $m2{$mID}{'cdsLines'}[-1][4];
        goto CHK_M_CDSLINES;
      }
    } elsif ($m2{$mID}{'cdsLines'}[0][6] eq '+' and $m2{$mID}{'cdsLines'}[0][7]  =~ m!^[^.0]$!) {
      scalar(@{$m2{$mID}{'5utrLines'}}) > 0 and die "[Err] FFFF 8Mp: [$mID] Should not have 5'-utr when the first codon is not intact.\n";
      $m2{$mID}{'exonLines'}[0][3] == $m2{$mID}{'cdsLines'}[0][3] or die "[Err] FFFF 7Mp: [$mID] additional 5'-exons.\n";
      $m2{$mID}{'cdsLines'}[0][3] += $m2{$mID}{'cdsLines'}[0][7];
      $m2{$mID}{'cdsLines'}[0][7] = 0;
      $m2{$mID}{'exonLines'}[0][3] = $m2{$mID}{'cdsLines'}[0][3];
      $m2{$mID}{'mrnaLine'}[3] = $m2{$mID}{'cdsLines'}[0][3];
      if ($m2{$mID}{'cdsLines'}[0][3] > $m2{$mID}{'cdsLines'}[0][4]) {
        # This CDS block should be removed entirely and the related exon should be removed too, because this exon should contain the CDS only.
        shift(@{$m2{$mID}{'cdsLines'}});
        shift(@{$m2{$mID}{'exonLines'}});
        $m2{$mID}{'mrnaLine'}[3] = $m2{$mID}{'cdsLines'}[0][3];
        goto CHK_M_CDSLINES;
      }
    }
  }# if (scalar(@{$m2{$mID}{'cdsLines'}}) > 0)
  my @o2;
  if ($m2{$mID}{'mrnaLine'}[6] eq '-') {
    @o2 = sort { $b->[4] <=> $a->[4] || $good_ele{lc($a->[2])} <=> $good_ele{lc($b->[2])} } (
      @{$m2{$mID}{'exonLines'}},@{$m2{$mID}{'5utrLines'}},@{$m2{$mID}{'cdsLines'}},@{$m2{$mID}{'3utrLines'}}
    );
  } else {
    @o2 = sort { $a->[3] <=> $b->[3] || $good_ele{lc($a->[2])} <=> $good_ele{lc($b->[2])} } (
      @{$m2{$mID}{'exonLines'}},@{$m2{$mID}{'5utrLines'}},@{$m2{$mID}{'cdsLines'}},@{$m2{$mID}{'3utrLines'}}
    );
  }
  print STDOUT join("\t", @{$m2{$mID}{'mrnaLine'}})."\n";
  for my $lR2 (@o2) {
    print STDOUT join("\t", @$lR2)."\n";
  }
  
}

sub load_gff3 {
  my ($fn) = @_;
  my %genes;
  my %mrnas;

  open(IGFF,'<',$fn) or die "FFFF 5\n";
  my $cntLineN = 0;
  while (<IGFF>) {
    $cntLineN ++;
    m!^\s*(#|$)! and do { next;};
    chomp;
    my @ta=split(/\t/, $_);
    my $eleID = lc($ta[2]);
    defined $good_ele{$eleID} or die "[Err] Don't recoganize element [$ta[2]]\n";
    my ($curID, $parentID);
    $ta[8] =~ m!(?:^|;)\s*ID\s*=\s*([^;\s]+)!i and $curID = $1;
    $ta[8] =~ m!(?:^|;)\s*Parent\s*=\s*([^;\s]+)! and $parentID = $1;
    if      ($eleID eq 'gene') {
      defined $genes{$curID} and defined $genes{$curID}{'geneLine'} and die "[Err] repeat geneID [$curID]: $_\n";
      $genes{$curID}{'geneLine'} = [@ta];
      $genes{$curID}{'lineN'} //= $cntLineN;
    } elsif ($eleID eq 'mrna') {
      defined $mrnas{$curID} and defined $mrnas{$curID}{'mrnaLine'} and die "[Err] Repeat mrnaID [$curID]: $_\n";
      $mrnas{$curID}{'mrnaLine'} = [@ta];
      $mrnas{$curID}{'geneID'} = $parentID;
      $genes{$parentID}{'mrnaIDs'}{$curID} = $cntLineN;
      $mrnas{$curID}{'lineN'} //= $cntLineN; # Don't want to use $.;
      $ta[6] =~ m!^[+-]$! or die "[Err] Bad strand in line: $_\n";
      $mrnas{$curID}{'str'} = $ta[6];
    } elsif ($eleID eq 'exon') {
      defined $parentID or die "[Err] No parentID found at line: $_\n";
      push(@{$mrnas{$parentID}{'exonLines'}}, [@ta]);
    } elsif ($eleID eq 'five_prime_utr') {
      defined $parentID or die "[Err] No parentID found at line: $_\n";
      push(@{$mrnas{$parentID}{'5utrLines'}}, [@ta]);
    } elsif ($eleID eq 'cds') {
      defined $parentID or die "[Err] No parentID found at line: $_\n";
      push(@{$mrnas{$parentID}{'cdsLines'}}, [@ta]);
    } elsif ($eleID eq 'three_prime_utr') {
      defined $parentID or die "[Err] No parentID found at line: $_\n";
      push(@{$mrnas{$parentID}{'3utrLines'}}, [@ta]);
    } else {
      die "[Err] Unknown element [$ta[2]]\n";
    }
  }
  close(IGFF);

  return(\%genes, \%mrnas);
}# load_gff3()

