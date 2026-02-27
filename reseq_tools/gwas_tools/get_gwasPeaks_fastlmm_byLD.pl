#!/usr/bin/perl
# 12/2/2024: Correct 'span' in .mergedblk file.
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 outPrefix minuslogP_cut in_fastlmm.res bo.vcf.gz peak_flank_distance_1e6 block_join_distance_100e3\n";

my $opref = shift;
my $pCut=shift;
my $fn=shift;
my $boVCF = '/data/Sunhh/wmhifi/analysis/sv_association/gwas_gbat/genotype/boAll.vcf.gz';
$boVCF = shift;
my $dist=1e6;
scalar(@ARGV) > 0 and $dist = shift;
my $sepDist = 1e5;
scalar(@ARGV) > 0 and $sepDist = shift; # If two blocks have distance <= this value, they are joined.

my $pl_peak2LD = 'perl /home/sunhh/tools/github/NGS_data_processing/reseq_tools/gwas_tools/peak2LD_plink_fastlmm.pl';

my $wdir = &fileSunhh::new_tmp_dir('create' => 1);

# Obtain P value cutoff if the input is a file (.cuts) instead of a number.
if (-f $pCut) {
  my @a1=&fileSunhh::load_tabFile($pCut);
  $a1[0][0] eq 'cut1_signif' or die "[Err] Bad input_pCut [$pCut]\n";
  $pCut = 10**(-$a1[1][0]);
} else {
  $pCut = 10**(-$pCut);
}

# Obtain markers with significant P values.
my $ifh=&openFH($fn);
my $ofh=&openFH("$fn.sig1", '>');
my @vars;
while (<$ifh>) {
  chomp;
  my @ta=split(/\t/,$_);
  # $ta[5] <= $pCut or next; # For emmax .p file.
  $ta[0] eq 'SNP' and next;
  $ta[4] <= $pCut or next; # For FaSTLMM 
  push(@vars, [@ta]); # In the format of FaSTLMM.res;
  print {$ofh} "$_\n";
}
close($ofh);
close($ifh);

# Define LD block for significant markers from top to bottom (sorted by P values);
my @blocks; # ([peakSNP, chr, start, end, span], [], ...); o1_LD.blk
open OPEAK,'>',"$opref.peak" or die;
open OBLK,'>',"$opref.blk" or die; # peakID chrV start end span;
print OBLK join("\t", qw/peakID chrV start end span/)."\n";
open OLINKVAR,'>',"$opref.tagvar" or die; # peakID chrV proxPosi r2 mrkID dist
print OLINKVAR join("\t", qw/peakID chrV proxPosi r2 mrkID dist/)."\n";
my %recorded_tagvar;

for my $l (sort {$a->[4] <=> $b->[4]} @vars) {
  my $is_inLD=0; # If this site has been included by a peak-based LD.
  for my $pt (@blocks) {
    if ($l->[1] eq $pt->[1] and $l->[3] >= $pt->[2] and $l->[3] <= $pt->[3]) {
      # Site is in peak block. $l->[1 & 3] in fastlmm, [0 & 1] in emmax result.
      $is_inLD = 1;
      # Check if this is recorded in tagVAR.
      if (!(defined $recorded_tagvar{$l->[0]})) {
        $recorded_tagvar{$l->[0]} = 1;
        print OLINKVAR join("\t", $pt->[0], $pt->[1], $l->[3], -1, $l->[0], 'NA')."\n";
      }
      last;
    }
  }
  if ($is_inLD == 0) {
    &fileSunhh::write2file("$wdir/i1", join("\t", @$l)."\n", ">");
    system "$pl_peak2LD $wdir/o1 $wdir/i1 $boVCF 0.8 $dist 1> /dev/null";
    # Record peaks and blocks.
    print OPEAK join("\t", @$l)."\n";
    open F1,'<',"$wdir/o1_LD.blk" or die;
    while (<F1>) {
      chomp;
      my @ta=split(/\t/, $_);
      $ta[0] eq 'peakID' and next;
      push(@blocks, [@ta]);
      print OBLK "$_\n";
    }
    close F1;
    # Record associated (linked) variants.
    open F2,'<',"$wdir/o1_LD.tagVAR" or die;
    while (<F2>) {
      chomp;
      my @ta=split(/\t/, $_);
      $ta[0] eq 'peakID' and next;
      print OLINKVAR "$_\n";
      $recorded_tagvar{$ta[4]} = 1;
    }
    close F2;
    system "rm -f $wdir/o1* $wdir/i1";
  }
}
close OPEAK;
close OBLK;
close OLINKVAR;

# Merge overlapping LD blocks
my @mergedblocks;
for (@blocks) { push(@mergedblocks, [@{$_}[0..4], ["[$_->[0]:$_->[1]:$_->[2]-$_->[3]:$_->[4]]"]]); }
my $prevMergedBlkN  = -1;
while ($prevMergedBlkN != scalar(@mergedblocks)) {
  $prevMergedBlkN = scalar(@mergedblocks);
  my @newmerged;
  for my $b1 (@mergedblocks) {
    my $is_ovl = 0;
    for my $b2 (@newmerged) {
      $b1->[1] eq $b2->[1] or next;
      $b1->[2] > $b2->[3]+$sepDist and next;
      $b1->[3] < $b2->[2]-$sepDist and next;
      $is_ovl = 1;
      push(@{$b2->[5]}, "[$b1->[0]:$b1->[1]:$b1->[2]-$b1->[3]:$b1->[4]]");
      $b1->[2] < $b2->[2] and $b2->[2] = $b1->[2];
      $b1->[3] > $b2->[3] and $b2->[3] = $b1->[3];
      $b2->[4] = $b2->[3]-$b2->[2]+1;
    }
    if ($is_ovl == 0) {
      push(@newmerged, [@$b1]);
    }
  }
  @mergedblocks = @newmerged;
}
open OMERGEBLK,'>',"$opref.mergedblk" or die;
print OMERGEBLK join("\t", qw/peakID chrV start end span containedPeakIDs/)."\n";
open OMAPBLKID,'>',"$opref.map_blkID" or die;
print OMAPBLKID join("\t", qw/old_peak new_peak/)."\n";
for (@mergedblocks) {
  print OMERGEBLK join("\t", @{$_}[0..4], join(";", @{$_->[5]}))."\n";
  for my $b3 (@{$_->[5]}) { $b3 =~ m!^\[(\S+):\S+:\d+\-\d+:\d+\]$! or die "ID=|$b3|\n"; print OMAPBLKID join("\t", $1, $_->[0])."\n"; }
}
close OMERGEBLK;
close OMAPBLKID;
system "deal_table.pl $opref.peak -kSrch_idx $opref.mergedblk > $opref.mergedpeak";
system "ColLink.pl $opref.tagvar -f1 $opref.map_blkID -Col2 0-5 -Col1 1 -fill merged_peakID | deal_table.pl -column 6,1-5,0 |deal_table.pl -best_uniqCol 4 -select_col 3 -select_rule 1|deal_table.pl -UniqColLine 4 > $opref.mergedtagvar";
system "rm -f $opref.tagvar $opref.peak";

&fileSunhh::_rmtree($wdir);

