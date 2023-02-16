#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;
use mathSunhh;
my $ms_obj = mathSunhh->new();

!@ARGV and die "perl $0  final.wiDropped.cds.blk.gff3  final_rmQloc.wiRepre.fmt.PAV.5cols  > final_rmQloc.wiRepre.fmt.PAV.rmOlap\n";

my $f1gff = shift;
my $f2tab = shift;

&tsmsg("[Log] Loading GFF\n");
my (%mrnaLoc);
%mrnaLoc = &loadGff3($f1gff);
#  $back{'byMID'} = \%loc;
#  $back{'byCHR'} = \%byChr;
&tsmsg("[Log] Loading groups\n");
my (%ogrp);
%{$ogrp{'grp2gene'}} = &loadGrp($f2tab); # {'grp2gene'}{$grpID} = [[g1_1,g1_2,g1_3], [g2_1,g2_2,..], [], []]
my %prefRank = qw(
GrpSyn_  1
GrpSep_  2
GrpOrp_  3
);
my @grpIDs = sort { $prefRank{substr($a,0,7)} <=> $prefRank{substr($b,0,7)} || $a cmp $b } keys %{$ogrp{'grp2gene'}};

# Get locations within each group.
&tsmsg("[Log] Getting locations\n");
for my $grpID (@grpIDs) {
  $ogrp{'grp2loc'}{$grpID} //= [];
  for (my $i=0;$i<@{$ogrp{'grp2gene'}{$grpID}};$i++) {
    $ogrp{'grp2loc'}{$grpID}[$i] = [];
    for (my $j=0;$j<@{$ogrp{'grp2gene'}{$grpID}[$i]};$j++) {
      my $mID = $ogrp{'grp2gene'}{$grpID}[$i][$j];
      defined $mrnaLoc{'byMID'}{$mID} or &stopErr("[Err] undefined mrnaID [$mID]\n");
      my $cdsSum = 0;
      for my $t1 (@{$mrnaLoc{'byMID'}{$mID}{'cds'}}) {
        $cdsSum += ($t1->[1]-$t1->[0]+1);
      }
      push(@{$ogrp{'grp2loc'}{$grpID}[$i]}, [$mrnaLoc{'byMID'}{$mID}{'chr'}, $mrnaLoc{'byMID'}{$mID}{'str'}, $mrnaLoc{'byMID'}{$mID}{'cds'}, $cdsSum]);
    }
  }
}

# Merge combined groups;
&tsmsg("[Log] Combining groups\n");
my %newGrp;
my %hasGrouped;
for (my $i=0;$i<@grpIDs;$i++) {
  $i % 10 == 0 and &tsmsg("[Log] Checking $i+1 th group [$grpIDs[$i]]\n");
  defined $hasGrouped{$grpIDs[$i]} and next;
  $hasGrouped{$grpIDs[$i]} = 1;
  my $refLocAH = $ogrp{'grp2loc'}{$grpIDs[$i]}; # [ _species 1, [chrID, str, [cds_se,...], cdsLen], [chrID_2, str2, [], cdsLen2], ... ]
  my $min_ovl_rate = 0.6;
  $newGrp{$grpIDs[$i]} = [$grpIDs[$i]];
  CHECK_GROUP1:
  for (my $j=$i+1;$j<@grpIDs;$j++) {
    defined $hasGrouped{$grpIDs[$j]} and next CHECK_GROUP1;
    my $chkLocAH = $ogrp{'grp2loc'}{$grpIDs[$j]};
    my $ovl_spec = 0;
    my $ttl_spec = 0;
    CHECK_SPECIES1:
    for (my $kR=0;$kR<@$refLocAH;$kR++) {
      # In each $kR species.
      scalar(@{$refLocAH->[$kR]}) == 0 and next;
      $ttl_spec ++;
      my $is_ovl = 0;
      CHECK_GENE1:
      for (my $m=0; $m<@{$refLocAH->[$kR]}; $m++) {
        # Compare each ref gene ($m) within same group and same specieus.
        for (my $n=0; $n<@{$chkLocAH->[$kR]}; $n++) {
          $refLocAH->[$kR][$m][0] eq $chkLocAH->[$kR][$n][0] or next;
          $refLocAH->[$kR][$m][1] eq $chkLocAH->[$kR][$n][1] or next;
          my ($ovlLen_nonDup, $ovlCnt_mayDup, $ovlLocR) = $ms_obj->compare_number_list($refLocAH->[$kR][$m][2], $chkLocAH->[$kR][$n][2], 'compare'=>'ovl', 'sort'=>0);
          # warn "$refLocAH->[$kR][$m][0]: $ovlLen_nonDup >= $min_ovl_rate * $refLocAH->[$kR][$m][3]\n";
          # warn "$chkLocAH->[$kR][$n][0]: $ovlLen_nonDup >= $min_ovl_rate * $chkLocAH->[$kR][$n][3]\n\n";
          if ($ovlLen_nonDup >= $min_ovl_rate * $refLocAH->[$kR][$m][3] 
            or $ovlLen_nonDup >= $min_ovl_rate * $chkLocAH->[$kR][$n][3]) {
            $is_ovl = 1;
            last;
          }
        }
        $is_ovl == 1 and last CHECK_GENE1;
      }# for check_gene1;
      $is_ovl == 1 and $ovl_spec ++;
    }# for check_species1;
    # warn "final ovl_spec=$ovl_spec\n";
    if ($ovl_spec >= $ttl_spec) {
      push(@{$newGrp{$grpIDs[$i]}}, $grpIDs[$j]);
      $hasGrouped{$grpIDs[$j]} = 1;
    }
  }# for check_group1
}

# Output groups;
my @grpIDs_new = sort { $prefRank{substr($a,0,7)} <=> $prefRank{substr($b,0,7)} || $a cmp $b } keys %newGrp;
for my $grpID (@grpIDs_new) {
  my @o1 = ($grpID);
  for (my $i=0;$i<@{$ogrp{'grp2gene'}{$grpID}};$i++) {
    my @ogenes;
    for my $id2 (@{$newGrp{$grpID}}) {
      push(@ogenes, @{$ogrp{'grp2gene'}{$id2}[$i]});
    }
    push(@o1, join(";", map { s!^\S\Span:!!; $_; } @ogenes));
  }
  print STDOUT join("\t", @o1)."\n";
}

#0       groupID GrpSyn_000001
#7       CLpan_geneIDs   Cla97C09G178432.1;Cla97C09G178434.1;Cla97C09G178436.1;Cla97C09G178438.1
#8       CMpan_geneIDs   CmUC09G176700.1;CmUC09G176710.1;CmUC09G176720.1;CmUC09G176730.1
#9       CApan_geneIDs   CaUC09G173540.1;CaUC09G173550.1;CaUC09G173640.1;CaUC09G173650.1;CaUC09G173660.1;CaUC09G173670.1
#10      CCpan_geneIDs   CcUC09G179820.1;CcUC09G179830.1;CcUC09G179840.1;CcUC09G179850.1


sub loadGrp {
  my ($f1) = @_;
  my (%loc);
  my @pref=qw(CLpan: CMpan: CApan: CCpan:);
  my $fh1 = &openFH($f1, '<');
  while (<$fh1>) {
    chomp;
    m!^\s*(#|$)! and next;
    my @ta=split(/\t/, $_);
    $ta[0] eq 'groupID' and next;
    for (my $i=0;$i<@pref;$i++) {
      my $j=$i+1;
      if (defined $ta[$j] and $ta[$j] !~ m!^\s*$!) {
        push(@{$loc{$ta[0]}}, [map { "$pref[$i]$_" } split(/;/, $ta[$j])]);
      } else {
        push(@{$loc{$ta[0]}}, []);
      }
    }
  }
  close($fh1);
  return (%loc);
}# loadGrp()

sub loadGff3 {
  my ($f1) = @_;
  my %loc;
  my %byChr;
  my $fh1 = &openFH($f1, '<');
  while (<$fh1>) {
    m!^\s*(#|$)! and next;
    chomp;
    my @ta=split(/\t/, $_);
    if ($ta[2] eq 'gene') {
      next;
    } elsif ($ta[2] eq 'mRNA') {
      next;
    } elsif ($ta[2] eq 'CDS') {
      $ta[8] =~ m!(?:^|\s)Parent=([^;]+)!i or &stopErr("[Err] unparsed ta8 [$ta[8]]\n");
      my $mID = $1;
      push(@{$loc{$mID}{'cds'}}, [@ta[3,4]]);
      $loc{$mID}{'chr'} //= $ta[0];
      $loc{$mID}{'chr'} eq $ta[0] or &stopErr("[Err] Diff CHR [$loc{$mID}{'chr'} VS. $ta[0]]\n");
      push(@{$byChr{$ta[0]}}, $mID);
      if ($ta[6] =~ m!^[+\-]$!) {
        $loc{$mID}{'str'} //= $ta[6];
        $loc{$mID}{'str'} eq $ta[6] or &stopErr("[Err] Why str ?[$ta[6]]\n");
      }
      $loc{$mID}{'se'}[0] //= $ta[3]; $loc{$mID}{'se'}[0] > $ta[3] and $loc{$mID}{'se'}[0] = $ta[3];
      $loc{$mID}{'se'}[1] //= $ta[4]; $loc{$mID}{'se'}[1] < $ta[4] and $loc{$mID}{'se'}[1] = $ta[4];
    } else {
      &stopErr("[Err] unparsed ele-type [$ta[2]]\n");
    }
  }
  close($fh1);
  for (keys %byChr) {
    @{$byChr{$_}} = sort { $loc{$a}{'se'}[0] <=> $loc{$b}{'se'}[0] || $loc{$a}{'se'}[1] <=> $loc{$b}{'se'}[1]  } @{$byChr{$_}};
  }
  for (keys %loc) {
    @{$loc{$_}{'cds'}} = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]} @{$loc{$_}{'cds'}};
  }
  my %back;
  $back{'byMID'} = \%loc;
  $back{'byCHR'} = \%byChr;
  return (%back);
}# loadGff3()


