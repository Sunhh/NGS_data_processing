#!/usr/bin/perl
# [5/25/2022]: Add help information. Given a list of gff3 files, retrieve the longest-CDS gene model at each non-overlaping region.
# [5/25/2022]: Add a rank to preferentially select the gene model in the gff file (which is the reference itself).
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;
use mathSunhh;
my $ms_obj = mathSunhh->new();

!@ARGV and die "[Err] perl $0 100 0.1 gff_list > slct_longestCDS.mrnaID_list\n";

my $pl_dealGff3 = 'perl /home/Sunhh/tools/github/NGS_data_processing/temp/deal_gff3.pl';
my $max_olapLen = shift; # 100
my $max_olapRat = shift; # 0.1

my $wrkDir = &fileSunhh::new_tmp_dir('create' => 1);
my @inGffFn;

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  push(@inGffFn, $ta[0]);
}

my %gene_models_onChr;
for (my $i=0; $i<@inGffFn; $i++) {
  &runCmd("$pl_dealGff3 -inGff $inGffFn[$i] -getJnLoc > $wrkDir/$i.gff3.JnLoc\n");
  my $rank = ($i == 0) ? 0 : 1;
  my $jnLoc_hr = &load_gff3JnLoc( $inGffFn[$i], "$wrkDir/$i.gff3.JnLoc", $rank );
  for my $chrID (keys %$jnLoc_hr) {
    for my $str (keys %{$jnLoc_hr->{$chrID}}) {
      push(@{$gene_models_onChr{$chrID}{$str}}, @{ $jnLoc_hr->{$chrID}{$str} }); # %$jnLoc_hr will not be released in the memory.
    }
  }
}

# Select the longest without overlapping.
for my $chrID (sort keys %gene_models_onChr) {
  for my $str (sort keys %{$gene_models_onChr{$chrID}}) {
    # Sort by CDS length;
    my @len_sort = sort { $b->[2] <=> $a->[2] || $a->[7] <=> $b->[7] || $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] || $a->[6] <=> $b->[6]} @{$gene_models_onChr{$chrID}{$str}};
    my @loc_sort = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] } @{$gene_models_onChr{$chrID}{$str}};
    # Find non-overlapping genes.
    my @slct_models;
    my %is_skip;
    for my $e_len (@len_sort) {
      # [ $tkey, $mrnaID, $LenInCDS, $CDSStart, $CDSEnd, \@cds_blocks, $CDSBlocksNum]
      defined $is_skip{$e_len->[0]} and next;
      push(@slct_models, $e_len);
      $is_skip{$e_len->[0]} = 1;
      my @new_loc_sort;
      for my $e_loc (@loc_sort) {
        defined $is_skip{$e_loc->[0]} and next;
        if      ($e_loc->[4] < $e_len->[3]) {
          push(@new_loc_sort, $e_loc);
        } elsif ($e_loc->[3] > $e_len->[4]) {
          push(@new_loc_sort, $e_loc);
        } else {
          my ($ovlLen_nonDup, $ovlCnt_mayDup, $ovlLoc_ar) = $ms_obj->compare_number_list($e_len->[5], $e_loc->[5], 'compare'=>'ovl', 'sort'=>0);
          if ($ovlLen_nonDup > $max_olapLen or $ovlLen_nonDup > $max_olapRat * $e_len->[2]) {
            # This gene ($e_loc) should be removed.
            $is_skip{$e_loc->[0]} = 1;
          } else {
            push(@new_loc_sort, $e_loc);
          }
        }
      }
      @loc_sort = @new_loc_sort;
    }
    for my $e_slct (@slct_models) {
      print STDOUT join("\t", $chrID, $str, @{$e_slct}[0,2,3,4])."\n";
    }
  }# for my $str (
}# for my $chrID (

&fileSunhh::_rmtree($wrkDir);

sub load_gff3JnLoc {
  my ($uniqID, $fn, $rank) = @_;
  $rank //= 0;
  my $back_hr = {};
  open F,'<',"$fn" or die;
  while (<F>) {
    chomp;
    my @ta=split(/\t/, $_);
    $ta[0] eq 'mrnaID' and next;
    my ($mrnaID, $SeqID, $mrnaStrand, $CDSStart, $CDSEnd, $LenInCDS, $CDSBlocks, $CDSBlocksNum) = @ta[0,2,5,6,7,8,9,10];
    my @cds_blocks;
    for my $t1 (split(/;/, $CDSBlocks)) {
      $t1 =~ m!^(\d+),(\d+)$! or &stopErr("[Err] Bad format of CDSBlocks [$t1]\n");
      my ($ts, $te) = ($1, $2);
      push(@cds_blocks, [$ts, $te]); # $ts is always <= $te, and the locations in @cds_blocks is ascending along elements disregarding with the strandness.
    }
    my $tkey = "$uniqID\t$mrnaID";
    push(@{$back_hr->{ $SeqID }{ $mrnaStrand }}, [ $tkey, $mrnaID, $LenInCDS, $CDSStart, $CDSEnd, \@cds_blocks, $CDSBlocksNum, $rank]); # {SeqID}{strand} = [ [], [], ... ]
  }
  close F;
  return($back_hr);
}
# header of a .gff3.JnLoc file.
# 0       mrnaID  CcU09G00640.1
# 1       geneID  CcU09G00640
# 2       SeqID   CmU531Chr09
# 3       mrnaStart       465096
# 4       mrnaEnd 469079
# 5       mrnaStrand      -
# 6       CDSStart        465407
# 7       CDSEnd  468915
# 8       LenInCDS        549
# 9       CDSBlocks       465407,465714;468502,468693;468867,468915
# 10      CDSBlocksNum    3
# 11      5UTR    164
# 12      3UTR    311
#
