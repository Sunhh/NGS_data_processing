#!/usr/bin/perl
# 5/26/2023: Updated to include 'tandem_duplication', 'simple relocation' (relocation), and 'collapsed_tandem_repeat'.
# 6/2/2023: Remove any VARs that are around 'reshuffling' boundaries (10 bp).
#           VARs around N gaps (10 bp) are also removed.
# 6/6/2023: I don't know why, but sometimes the reference base is same as query base in substitution.
# 6/8/2023: Shift <INV> REF position to 1 bp left in VCF.
# 8/15/2023: Accept SVs around 'reshuffling' segments, because I found true relocation-type DELs aside reshufflings.
#            This can introduce many false positives in INS!!!
# 8/16/2023: Return to previous version at 6/8/2023.
# 8/16/2023: Change 'SVTYPE' information to ease further check.
# 8/16/2023: Try to include query location for further check.
# 8/17/2023: Trim substitutions with same leading/tailing bases.
# 9/8/20923: Fixing a bug for 'relocation-insertion' on reverse strand.

use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 ref.fa.gz qry.fa.gz sampleID out_ref_struct.gff.gz > out_ref_struct.vcf\n";

my $fnRFa = shift;
my $fnQFa = shift;
my $samID = shift;
my $fnGff = shift;

my $maxDist = 10;
my $wdir = &fileSunhh::new_tmp_dir('create' => 1);
my $pl_fmtPaf = "/home/Sunhh/tools/github/NGS_data_processing/evolution_tools/SV_detection/fmt_paf.pl";
my $exe_stretcher = "stretcher";


# Possible Names ($svType):
###  Accept INS+DEL: substitution;
###  Accept INS: insertion/duplication/tandem_duplication;
###  Accept INS: relocation-insertion;
###  Accept DEL: deletion/collapsed_repeat/collapsed_tandem_repeat;
###  Accept DEL: relocation;
###  Accept INV: inversion;
###  Ignore: gap/inserted_gap/unaligned_end/unaligned_beginning/...;
###  Ignore: relocation*/translocation*/inversion/reshuffling/...;

#
# Basic settings.
#
my %goodName;
for (qw/deletion collapsed_repeat insertion duplication relocation-insertion inversion substitution/) {$goodName{$_}=1;}
$goodName{'tandem_duplication'} = 1;
$goodName{'collapsed_tandem_repeat'} = 1;
$goodName{'relocation'} = 1;
my %badName; # reshuffuling should be removed too.
for (qw/translocation translocation-insertion translocation-inserted_gap translocation-insertion_ATGCN translocation-overlap/) {$badName{$_} = 1;}
for (qw/relocation-inserted_gap relocation-insertion_ATGCN relocation-overlap/) {$badName{$_} = 1;}
for (qw/unaligned_beginning unaligned_end gap inserted_gap/) {$badName{$_} = 1;}
my %otherName;
my %restName;
my %bad_loc;

my %refSeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnRFa)};
for (keys %refSeq) { $refSeq{$_}{'seq'} =~ s!\s!!g; $refSeq{$_}{'len'}=length($refSeq{$_}{'seq'}); $refSeq{$_}{'seq'} = uc($refSeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded ref fa [$fnRFa]\n");

my %qrySeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnQFa)};
for (keys %qrySeq) { $qrySeq{$_}{'seq'} =~ s!\s!!g; $qrySeq{$_}{'len'}=length($qrySeq{$_}{'seq'}); $qrySeq{$_}{'seq'} = uc($qrySeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded qry fa [$fnQFa]\n");

{
# Add N gap list to %bad_loc;
my $ntxt = '[nNuU]+';
for my $k1 (keys %refSeq) {
  for my $a1 (&siteList(\$ntxt, \$refSeq{$k1}{'seq'}, 'Min')) {
    push(@{$bad_loc{'r'}{$k1}}, [$a1->[0], $a1->[1]]);
  }
}
for my $k1 (keys %qrySeq) {
  for my $a1 (&siteList(\$ntxt, \$qrySeq{$k1}{'seq'}, 'Min')) {
    push(@{$bad_loc{'q'}{$k1}}, [$a1->[0], $a1->[1]]);
  }
}
}

my $fhGff = &openFH($fnGff);
print STDOUT <<'HHHS';
##fileformat=VCFv4.2
##ALT=<ID=BND,Description="Translocation, which should be reshuffling in NucDiff.">
##ALT=<ID=CNV,Description="Copy Number Polymorphism">
##ALT=<ID=CPX,Description="Complex SV">
##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DEL:REP,Description="Deletion and collapsed repeat">
##ALT=<ID=DEL:TAN,Descirption="Collapsed tandem repeat">
##ALT=<ID=DEL:REL,Description="Relocation-shaped deletion">
##ALT=<ID=DEL:SUB,Description="Deletion-type substitution">
##ALT=<ID=DEL:ALN,Description="Deletion-type substitution from stretcher alignment.">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INS:DUP,Description="Duplication and/or insertion">
##ALT=<ID=INS:TAN,Description="Tandem duplication">
##ALT=<ID=INS:REL,Description="Relocation-insertion">
##ALT=<ID=INS:SUB,Description="Insertion-type substitution">
##ALT=<ID=INS:ALN,Description="Deletion-type substitution from stretcher alignment.">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">
##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">
##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">
##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">
##ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">
##ALT=<ID=UNK:SUB,Description="Mostly equal length substitutions from ref_struct">
##ALT=<ID=UNK:ALN,Description="Mostly equal length substitutions from re-aligned ref_struct">
##ALT=<ID=INV,Description="Inversion">
##CPX_TYPE_NDstruct_DEL="Deletion in ref_struct.gff generated from NucDiff"
##CPX_TYPE_NDstruct_INS="Insertion in ref_struct.gff generated from NucDiff"
##CPX_TYPE_NDstruct_SUB="Substitution in ref_struct.gff generated from NucDiff"
##CPX_TYPE_NDstruct_uALN="unaligned_beginning/unaligned_end region in ref_struct.gff generated from NucDiff"
##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Class of complex variant.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=ALTPOS,Number=1,Type=String,Description="Corresponding ALT sequence position in ChrID:start-end format. This is experimental">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
HHHS

for (sort keys %refSeq) {
  print STDOUT "##contig=<ID=$_,length=$refSeq{$_}{'len'}>\n";
}

print STDOUT join("\t", "#CHROM", qw/POS ID REF ALT QUAL FILTER INFO FORMAT/, $samID)."\n";
my %svPool; # {svType}{refID}{svID} = [refS, refE, qryID, qryS, qryE, strand];
my %allRefID;
while (<$fhGff>) {
  chomp;
  m!^\s*#! and next;
  my @ta=split(/\t/, $_);
  my @inf;
  push(@inf, "ALGORITHMS=NucDiff_struct");
  my %h8 = %{ &sepTA8($ta[8]) };
  defined $badName{$h8{'Name'}} and next;
  defined $otherName{$h8{'Name'}} and do { &tsmsg("[Wrn] Skip Name [$h8{'Name'}].\n"); next; };
  if ($h8{'Name'} =~ m!^reshuffling!) {
    push(@{$bad_loc{'r'}{$ta[0]}}, [@ta[3,4]]);
    $h8{'query_coord'} =~ m!^(\d+)\-(\d+)$! or die "bad 1: $_\n";
    my ($qs,$qe) = ($1, $2);
    push(@{$bad_loc{'q'}{$h8{'query_sequence'}}}, [$qs, $qe]);
    next;
  }
  defined $goodName{$h8{'Name'}} or do { &tsmsg("[Wrn] Skip unrecognized Name [$h8{'Name'}]\n"); next; };

  my ($refID, $refS, $refE);
  my ($qryID, $qryS, $qryE);
  my ($strN, $svType, $svID) = (1, $h8{'Name'}, $h8{'ID'});
  $refID = $ta[0];
  $qryID = $h8{'query_sequence'};
  if ($svType eq 'deletion') {
    $refS=$ta[3]; $refE=$ta[4]; $strN = $h8{'query_dir'};
    $qryS = $qryE = $h8{'query_coord'}; # When strN=-1, the inverted_ref is inserted on the right of qryS.
  } elsif ($svType eq 'collapsed_repeat') {
    unless (defined $h8{'ref_repeated_region'} and $h8{'ref_repeated_region'} ne '') {
      # When the following SV is 'reshuffling', there is no 'ref_repeated_region' tag.
      next;
    }
    $strN = $h8{'query_dir'};
    $qryS = $qryE = $h8{'query_coord'}; # dir=-1 doesn't change qryS/qryE;
    $h8{'ref_repeated_region'} =~ m!^(\d+)\-(\d+)$! or die "Fk 2:$_\n";
    $refS = $2+1; $refE = $ta[4];
    if ($strN == -1) { $refS = $ta[3]; $refE = $1-1; }
  } elsif ($svType eq 'insertion') {
    $refS=$refE=$ta[3]; $strN = $h8{'query_dir'};
    $h8{'query_coord'} =~ m!^(\d+)\-(\d+)$! or die "Fk 4\n";
    $qryS = $1; $qryE = $2; # dir=-1 doesn't change refS.
  } elsif ($svType eq 'duplication') {
    $refS=$refE=$ta[3]; $strN = $h8{'query_dir'}; # refS/refE is not affected by strand.
    $h8{'query_coord'} =~ m!^(\d+)\-(\d+)$! or die "Fk 5\n";
    $qryS = $1; $qryE = $2; # !!! $qryS needs to be updated in neighboring insertion SV no matter on what strand it aligned.
  } elsif ($svType eq 'relocation-insertion') {
    # Currently the strand is not recorded. I should compare SV order or some other methods to determine it.
    $strN = 0;
    $h8{'blk_ref'} =~ m!^(\d+)\-(\d+)$! or die "Fk 3\n";
    ($refS, $refE) = ($1, $2); # They should be updated later, which should be always the left flank's end position.
    $h8{'breakpoint_query'} =~ m!^(\d+)\-(\d+)$! or die "Fk 6\n";
    ($qryS, $qryE) = ($1, $2);
  } elsif ($svType eq 'relocation') {
    $strN = 0;
    $h8{'blk_ref'} =~ m!^(\d+)\-(\d+)$! or die "Fk 3a\n";
    ($refS, $refE) = ($1, $2); # They should be updated later, which should be always the left flank's end position.
    $h8{'breakpoint_query'} =~ m!^(\d+)\-(\d+)$! or die "Fk 6a\n";
    ($qryS, $qryE) = ($1, $2);
  } elsif ($svType eq 'inversion') {
    $refS=$ta[3]; $refE=$ta[4]; $strN=$h8{'query_dir'}; # The strand doesn't affect positions.
    $h8{'query_coord'} =~ m!^(\d+)\-(\d+)$! or die "Fk 7\n";
    ($qryS, $qryE) = ($1, $2);
  } elsif ($svType eq 'substitution') {
    ($refS, $refE) = @ta[3,4]; $strN=$h8{'query_dir'};
    $h8{'query_coord'} =~ m!^(\d+)\-(\d+)$! or die "Fk 8\n";
    ($qryS, $qryE) = ($1, $2);
  } elsif ($svType eq 'tandem_duplication') {
    $refS=$refE=$ta[3]; $strN=$h8{'query_dir'};
    $h8{'query_coord'} =~ m!^(\d+)\-(\d+)$! or die "Fk 10\n"; # strand should not matter.
    ($qryS, $qryE) = ($1, $2);
  } elsif ($svType eq 'collapsed_tandem_repeat') {
    $refS=$ta[3]; $refE=$ta[4]; $strN=$h8{'query_dir'};
    $qryS=$qryE=$h8{'query_coord'};
  } else {
    defined $restName{$h8{'Name'}} or &tsmsg("[Wrn] Skip Name [$h8{'Name'}].\n");
    next;
  }
  defined $svPool{$svType}{$refID}{$svID} and die "Fk 9\n";
  $svPool{$svType}{$refID}{$svID} = [$refS, $refE, $qryID, $qryS, $qryE, $strN];
  $allRefID{$refID} //= 1;
}
close ($fhGff);
# Sort svID by position.
my %refID2svID; # refID2svID{svType}{refID} = [ svID_1, svID_2, ... ];
for my $svType (keys %svPool) {
  for my $refID (keys %{$svPool{$svType}}) {
    @{$refID2svID{$svType}{$refID}} = sort { $svPool{$svType}{$refID}{$a}[0] <=> $svPool{$svType}{$refID}{$b}[0] 
        || $svPool{$svType}{$refID}{$a}[1] <=> $svPool{$svType}{$refID}{$b}[1]
      } keys %{$svPool{$svType}{$refID}};
  }
}

#
# Deal with special cases in which two consecutive records are used together.
#
# Type combinations:
#   DEL:REP : deletion     + collapsed_repeat
#   DEL:REL : relocation: SV_XX.1 + SV_XX.2.
#   INS:DUP : insertion    + duplication
#   DEL:SUB : substitution + deletion  : Re-align.
#   INS:SUB : substitution + insertion : Re-align.
#   INS:REL : 'relocation-insertion': SV_XX.1 + SV_XX.2
my %new_svType;
for my $refID (sort keys %allRefID) {
  my (@del_svIDs, @cor_svIDs, @sub_svIDs, @dup_svIDs, @ins_svIDs, @rin_svIDs, @rel_svIDs);
  defined $refID2svID{'deletion'}{$refID} and @del_svIDs = @{$refID2svID{'deletion'}{$refID}};
  defined $refID2svID{'collapsed_repeat'}{$refID} and @cor_svIDs = @{$refID2svID{'collapsed_repeat'}{$refID}};
  defined $refID2svID{'substitution'}{$refID} and @sub_svIDs = @{$refID2svID{'substitution'}{$refID}};
  defined $refID2svID{'duplication'}{$refID} and @dup_svIDs = @{$refID2svID{'duplication'}{$refID}};
  defined $refID2svID{'insertion'}{$refID} and @ins_svIDs = @{$refID2svID{'insertion'}{$refID}};
  defined $refID2svID{'relocation-insertion'}{$refID} and @rin_svIDs = grep { $_ =~ m!\.1$! } @{$refID2svID{'relocation-insertion'}{$refID}};
  @rin_svIDs = map { s!\.1!!; $_; } @rin_svIDs;
  defined $refID2svID{'relocation'}{$refID} and @rel_svIDs = grep { $_ =~ m!\.1$! } @{$refID2svID{'relocation'}{$refID}};
  @rel_svIDs = map { s!\.1!!; $_; } @rel_svIDs;
  my %rm_del_idx;
  my %rm_ins_idx;
  # DEL:REP : deletion     + collapsed_repeat
  for my $cor_ID (@cor_svIDs) {
    my ($cor_rS, $cor_rE, $cor_qID, $cor_qS, $cor_qE, $cor_str) = @{$svPool{'collapsed_repeat'}{$refID}{$cor_ID}};
    for (my $i=0; $i<@del_svIDs; $i++) {
      defined $rm_del_idx{$i} and next;
      my $del_ID = $del_svIDs[$i];
      my ($del_rS, $del_rE, $del_qID, $del_qS, $del_qE, $del_str) = @{$svPool{'deletion'}{$refID}{$del_ID}};
      $del_rE < $cor_rS and next;
      $del_rS > $cor_rE and last;
      $cor_str == 1  and $del_str == 1  and $del_rS == $cor_rS and do {$rm_del_idx{$i}=1; $new_svType{$cor_ID}='DEL:REP'; last;};
      $cor_str == -1 and $del_str == -1 and $del_rE == $cor_rE and do {$rm_del_idx{$i}=1; $new_svType{$cor_ID}='DEL:REP'; last;};
    }
  }
  # DEL:REL : 'relocation': SV_XX.1 + SV_XX.2. Update location and strand.
  for my $rel_ID (@rel_svIDs) {
    my $rel_ID1 = "$rel_ID.1";
    my $rel_ID2 = "$rel_ID.2";
    my ($r1_rS, $r1_rE, $r1_qID, $r1_qS, $r1_qE, $r1_str) = @{$svPool{'relocation'}{$refID}{$rel_ID1}};
    my ($r2_rS, $r2_rE, $r2_qID, $r2_qS, $r2_qE, $r2_str) = @{$svPool{'relocation'}{$refID}{$rel_ID2}};
    $new_svType{$rel_ID1} = 'DEL:REL';
    $new_svType{$rel_ID2} = 'DEL:REL';
    if ($r1_rS > $r2_rS) {
      # Deleted ref region should be [$r2_rE+1, $r1_rS-1]
      $svPool{'relocation'}{$refID}{$rel_ID1}[5] = -1;
      $svPool{'relocation'}{$refID}{$rel_ID2}[5] = -1;
    } else {
      # Deleted ref region should be [$r1_rE+1, $r2_rS-1]
      $svPool{'relocation'}{$refID}{$rel_ID1}[5] = 1;
      $svPool{'relocation'}{$refID}{$rel_ID2}[5] = 1;
    }
  }

  
  # INS:DUP : insertion    + duplication
  for my $dup_ID (@dup_svIDs) {
    my ($dup_rS, $dup_rE, $dup_qID, $dup_qS, $dup_qE, $dup_str) = @{$svPool{'duplication'}{$refID}{$dup_ID}};
    for (my $i=0; $i<@ins_svIDs; $i++) {
      defined $rm_ins_idx{$i} and next;
      my $ins_ID = $ins_svIDs[$i];
      my ($ins_rS, $ins_rE, $ins_qID, $ins_qS, $ins_qE, $ins_str) = @{$svPool{'insertion'}{$refID}{$ins_ID}};
      $ins_rE < $dup_rS and next;
      $ins_rS > $dup_rE and last;
      if ( $ins_rS == $dup_rS and $dup_str*$ins_str == 1 ) {
        $rm_ins_idx{$i}=1;
        $svPool{'duplication'}{$refID}{$dup_ID}[3]=$ins_qS;
        $new_svType{$dup_ID} = 'INS:DUP';
        last;
      }
    }
  }

  CHK_COMB_SUB:
  for my $sub_ID (@sub_svIDs) {
    my ($sub_rS, $sub_rE, $sub_qID, $sub_qS, $sub_qE, $sub_str) = @{$svPool{'substitution'}{$refID}{$sub_ID}};
    # DEL:SUB : substitution + deletion  : Re-align. OR record directly. Re-align them later.
    for (my $i=0; $i<@del_svIDs; $i++) {
      defined $rm_del_idx{$i} and next;
      my $del_ID = $del_svIDs[$i];
      my ($del_rS, $del_rE, $del_qID, $del_qS, $del_qE, $del_str) = @{$svPool{'deletion'}{$refID}{$del_ID}};
      $del_rS > $sub_rE+1 and last;
      if ($sub_str == 1 and $del_str == 1 and $del_rS == $sub_rE+1) {
        $rm_del_idx{$i}=1;
        $svPool{'substitution'}{$refID}{$sub_ID}[1] = $del_rE;
        $new_svType{$sub_ID} = 'DEL:SUB';
        next CHK_COMB_SUB;
      } elsif ($sub_str == -1 and $del_str == -1 and $del_rE == $sub_rS-1) {
        $rm_del_idx{$i}=1;
        $svPool{'substitution'}{$refID}{$sub_ID}[0] = $del_rS;
        $new_svType{$sub_ID} = 'DEL:SUB';
        next CHK_COMB_SUB;
      }
    }
    # INS:SUB : substitution + insertion : Re-align. OR record long sequence. Re-align them later in VCF files.
    for (my $i=0; $i<@ins_svIDs; $i++) {
      defined $rm_ins_idx{$i} and next;
      my $ins_ID = $ins_svIDs[$i];
      my ($ins_rS, $ins_rE, $ins_qID, $ins_qS, $ins_qE, $ins_str) = @{$svPool{'insertion'}{$refID}{$ins_ID}};
      $ins_rS > $sub_rE and last;
      if ($sub_str == 1 and $ins_str == 1 and $ins_rS == $sub_rE) {
        $rm_ins_idx{$i}=1;
        $svPool{'substitution'}{$refID}{$sub_ID}[4] = $ins_qE;
        $new_svType{$sub_ID} = 'INS:SUB';
        next CHK_COMB_SUB;
      } elsif ($sub_str == -1 and $ins_str == -1 and $ins_rE+1 == $sub_rS) {
        $rm_ins_idx{$i}=1;
        $svPool{'substitution'}{$refID}{$sub_ID}[4] = $ins_qE;
        $new_svType{$sub_ID} = 'INS:SUB';
        next CHK_COMB_SUB;
      }
    }
  }
  # INS:REL : 'relocation-insertion': SV_XX.1 + SV_XX.2. Update location and strand.
  for my $rin_ID (@rin_svIDs) {
    my $rin_ID1 = "$rin_ID.1";
    my $rin_ID2 = "$rin_ID.2";
    my ($r1_rS, $r1_rE, $r1_qID, $r1_qS, $r1_qE, $r1_str) = @{$svPool{'relocation-insertion'}{$refID}{$rin_ID1}};
    my ($r2_rS, $r2_rE, $r2_qID, $r2_qS, $r2_qE, $r2_str) = @{$svPool{'relocation-insertion'}{$refID}{$rin_ID2}};
    $new_svType{$rin_ID1} = 'INS:REL';
    $new_svType{$rin_ID2} = 'INS:REL';
    if ($r1_rS > $r2_rS) {
      $svPool{'relocation-insertion'}{$refID}{$rin_ID1}[5] = -1;
      $svPool{'relocation-insertion'}{$refID}{$rin_ID2}[5] = -1;
    } else {
      $svPool{'relocation-insertion'}{$refID}{$rin_ID1}[5] = 1;
      $svPool{'relocation-insertion'}{$refID}{$rin_ID2}[5] = 1;
    }
  }
  # Clean deletion SVs.
  my (@new_del, @new_ins);
  for (my $i=0; $i<@del_svIDs; $i++) {
    if (defined $rm_del_idx{$i}) {
      delete $svPool{'deletion'}{$refID}{$del_svIDs[$i]};
    } else {
      push(@new_del, $del_svIDs[$i]);
    }
  }
  @{$refID2svID{'deletion'}{$refID}} = @new_del;
  # Clean insertion SVs.
  for (my $i=0; $i<@ins_svIDs; $i++) {
    if (defined $rm_ins_idx{$i}) {
      delete $svPool{'insertion'}{$refID}{$ins_svIDs[$i]};
    } else {
      push(@new_ins, $ins_svIDs[$i]);
    }
  }
  @{$refID2svID{'insertion'}{$refID}} = @new_ins;
}

#
# Process all VARs and output the VCF file.
#
for my $refID (sort keys %allRefID) {
  # Output deletion
  for my $var_ID (@{$refID2svID{'deletion'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'deletion'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-2, $rE-$rS+2);
    my $baseAlt = substr($baseRef, 0, 1);
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=DEL;END=$rE;SVLEN=-".($rE-$rS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS-1, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output collapsed_repeat (DEL:REP)
  for my $var_ID (@{$refID2svID{'collapsed_repeat'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'collapsed_repeat'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-2, $rE-$rS+2);
    my $baseAlt = substr($baseRef, 0, 1);
    my $t1 = (defined $new_svType{$var_ID}) ? $new_svType{$var_ID} : 'DEL:REP:UNK' ;
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=$t1;END=$rE;SVLEN=-".($rE-$rS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS-1, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output collapsed_tandem_repeat (DEL:TAN)
  for my $var_ID (@{$refID2svID{'collapsed_tandem_repeat'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'collapsed_tandem_repeat'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-2, $rE-$rS+2);
    my $baseAlt = substr($baseRef, 0, 1);
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=DEL:TAN;END=$rE;SVLEN=-".($rE-$rS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS-1, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output relocation (DEL:REL)
  for my $var_ID (grep { $_ =~ m!\.1$! } @{$refID2svID{'relocation'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'relocation'}{$refID}{$var_ID}};
    my $rel_ID2 = $var_ID; $rel_ID2 =~ s!\.1$!.2!;
    my ($rS2, $rE2, $qryID2, $qS2, $qE2, $strN2) = @{$svPool{'relocation'}{$refID}{$rel_ID2}};
    my ($rSx, $rEx) = ($rE+1, $rS2-1);
    $strN == -1 and ($rSx, $rEx) = ($rE2+1, $rS-1);
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rSx-2, $rEx-$rSx+2);
    my $baseAlt = substr($baseRef, 0, 1);
    my $t1 = (defined $new_svType{$var_ID}) ? $new_svType{$var_ID} : 'DEL:REL:UNK' ;
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=$t1;END=$rEx;SVLEN=-".($rEx-$rSx+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rSx-1, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output insertion (INS)
  for my $var_ID (@{$refID2svID{'insertion'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'insertion'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
    my $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
    $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
    $baseAlt = $baseRef . $baseAlt;
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=INS;END=".($rS+1).";SVLEN=".($qE-$qS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output duplication (INS:DUP)
  for my $var_ID (@{$refID2svID{'duplication'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'duplication'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
    my $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
    $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
    $baseAlt = $baseRef . $baseAlt;
    my $t1 = (defined $new_svType{$var_ID}) ? $new_svType{$var_ID} : 'INS:DUP:UNK' ;
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=$t1;END=".($rS+1).";SVLEN=".($qE-$qS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output tandem_duplication (INS:TAN)
  for my $var_ID (@{$refID2svID{'tandem_duplication'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'tandem_duplication'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
    my $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
    $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
    $baseAlt = $baseRef . $baseAlt;
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=INS:TAN;END=".($rS+1).";SVLEN=".($qE-$qS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output relocation-insertion (INS:REL)
  for my $var_ID (grep { $_ =~ m!\.1$! } @{$refID2svID{'relocation-insertion'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'relocation-insertion'}{$refID}{$var_ID}};
    my ($rS1, $rE1, $qryID1, $qS1, $qE1, $strN1) = ($rS, $rE, $qryID, $qS, $qE, $strN);
    my $rin_ID2 = $var_ID; $rin_ID2 =~ s!\.1$!.2!;
    my ($rS2, $rE2, $qryID2, $qS2, $qE2, $strN2) = @{$svPool{'relocation-insertion'}{$refID}{$rin_ID2}};
    $strN == -1 and ($rS, $rE, $qryID, $qS, $qE, $strN) = ($rE2, $rE2, $qryID2, $qS2, $qE2, $strN2);
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
    my $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
    $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
    $baseAlt = $baseRef . $baseAlt;
    my $t1 = (defined $new_svType{$var_ID}) ? $new_svType{$var_ID} : 'INS:REL:UNK' ;
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=$t1;END=".($rS+1).";SVLEN=".($qE-$qS+1).";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS1, $rE1, $qryID, $qS1, $qE1, $maxDist) and next;
    &boundary_inBad($refID, $rS2, $rE2, $qryID, $qS2, $qE2, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output inversion (INV)
  for my $var_ID (@{$refID2svID{'inversion'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'inversion'}{$refID}{$var_ID}};
    # my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, $rE-$rS+1);
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-2, 1); # Inverted region is [$rS, $rE], labeled at $rS-1.
    my $baseAlt = '<INV>';
    my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=INV;END=$rE;SVLEN=0".";ALTPOS=$qryID:$qS-$qE";
    &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
    print STDOUT join("\t", $refID, $rS-1, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
  }
  # Output substitution (DEL:SUB / INS:SUB)
  for my $var_ID (@{$refID2svID{'substitution'}{$refID}}) {
    my ($rS, $rE, $qryID, $qS, $qE, $strN) = @{$svPool{'substitution'}{$refID}{$var_ID}};
    my $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, $rE-$rS+1);
    my $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
    $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
    $baseRef eq $baseAlt and next;
    my $t1 = (defined $new_svType{$var_ID}) ? $new_svType{$var_ID} : 'UNK:SUB' ;
    # Try if minimap2 can map the them. If yes, go EMBOSS-stretcher. If no, keep it and remove leading/tailing same bases.
    &fileSunhh::write2file("$wdir/r.fa", ">r\n$baseRef\n", '>');
    &fileSunhh::write2file("$wdir/q.fa", ">q\n$baseAlt\n", '>');
    my $is_good = 0;
    open F1, '-|', "minimap2 $wdir/r.fa $wdir/q.fa -c 2> $wdir/stderr | perl $pl_fmtPaf" or &stopErr("[Err] Failed CMD: $!\n");
    while (my $l1 = <F1>) {
      chomp($l1);
      my @tb=split(/\t/, $l1);
      $tb[0] eq 'QID' and next;
      ($tb[12] >= 80 and $tb[3]-$tb[2] >= 0.6 * $tb[1]) or 
        ($tb[13] >= 80 and $tb[8]-$tb[7] >= 0.6 * $tb[6]) or next;
      $is_good = 1;
      last;
    }
    close F1;
    if ($is_good == 0) {
      my $pi=0;
      while ( substr($baseRef, $pi, 1) eq substr($baseAlt, $pi, 1) and $pi < length($baseRef) and $pi < length($baseAlt)) {
        $pi ++; $rS ++;
        if ($strN == -1) {
          $qE --;
        } else {
          $qS ++;
        }
      }
      if ($rE+1 == $rS) {
        ($rS, $rE) = ($rE, $rS);
        $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
        $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
        $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
        $baseAlt = $baseRef . $baseAlt;
      } elsif ($qE+1 == $qS) {
        ($qS, $qE) = ($qE, $qS);
        $rS --; 
        $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, $rE-$rS+1);
        $baseAlt = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
      } elsif ($rE >= $rS and $qE >= $qS) {
        $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, $rE-$rS+1);
        $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
        $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
        $pi = -1;
        while ( substr($baseRef, $pi, 1) eq substr($baseAlt, $pi, 1) and $pi*-1 > length($baseRef) and $pi*-1 > length($baseAlt)) {
          $pi--; $rE--;
          if ($strN == -1) {
            $qS ++;
          } else {
            $qE --;
          }
        }
        if ($rE+1 == $rS) {
          ($rS, $rE) = ($rE, $rS);
          $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
          $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
          $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
          $baseAlt = $baseRef . $baseAlt;
        } elsif ($qE+1 == $qS) {
          ($qS, $qE) = ($qE, $qS);
          $rS --; 
          $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, $rE-$rS+1);
          $baseAlt = substr($refSeq{$refID}{'seq'}, $rS-1, 1);
        } elsif ($rE >= $rS and $qE >= $qS) {
          $baseRef = substr($refSeq{$refID}{'seq'}, $rS-1, $rE-$rS+1);
          $baseAlt = substr($qrySeq{$qryID}{'seq'}, $qS-1, $qE-$qS+1);
          $strN == -1 and &fastaSunhh::rcSeq(\$baseAlt, 'rc');
        } else {
          &stopErr("[Err] Should not reach here 2: $var_ID $rS $rE $qS $qE\n");
        }
      } else {
        &stopErr("[Err] Should not reach here 1: $var_ID $rS $rE $qS $qE\n");
      }
      my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=$t1;END=$rE;SVLEN=".($qE-$qS-$rE+$rS).";ALTPOS=$qryID:$qS-$qE";
      &boundary_inBad($refID, $rS, $rE, $qryID, $qS, $qE, $maxDist) and next;
      print STDOUT join("\t", $refID, $rS, '.', $baseRef, $baseAlt, '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
    } else {
      # Use EMBOSS-stretcher to re-calculate variants.
      system("$exe_stretcher -asequence $wdir/r.fa -bsequence $wdir/q.fa -aformat fasta -outfile $wdir/stretcher.fa");
      my %alnSeq = %{ $fs_obj->save_seq_to_hash('faFile' => "$wdir/stretcher.fa") };
      $alnSeq{'r'}{'seq'} =~ s!\s!!g; my @rseq = split(//, uc($alnSeq{'r'}{'seq'}));
      $alnSeq{'q'}{'seq'} =~ s!\s!!g; my @qseq = split(//, uc($alnSeq{'q'}{'seq'}));
      my @prev = ('', '', '', 0, 0, 0, 0); # ['type', 'rseq', 'qseq', rS, rE, qS, qE] 1-based position;
      my @blks;
      for (my $pi=0; $pi<scalar(@rseq); $pi++) {
        if ($rseq[$pi] eq $qseq[$pi]) {
          if ($prev[0] eq 'same') {
            $prev[4] ++; $prev[6] ++;
          } else {
            push(@blks, [@prev]);
            @prev = ('same', '', '', $prev[4]+1, $prev[4]+1, $prev[6]+1, $prev[6]+1);
          }
        } elsif ($rseq[$pi] eq '-') {
          if ($prev[0] eq 'ins') {
            # Insert on the right of labeled position (rS==rE).
            $prev[2] .= $qseq[$pi]; $prev[6] ++;
          } else {
            push(@blks, [@prev]);
            @prev = ('ins', '', $qseq[$pi], $prev[4], $prev[4], $prev[6]+1, $prev[6]+1);
          }
        } elsif ($qseq[$pi] eq '-') {
          if ($prev[0] eq 'del') {
            # The QPOS should be always on the left of deleted RSEQ.
            $prev[1] .= $rseq[$pi]; $prev[4] ++;
          } else {
            push(@blks, [@prev]);
            @prev = ('del', $rseq[$pi], '', $prev[4]+1, $prev[4]+1, $prev[6], $prev[6]);
          }
        } else {
          # Substitution (may have Ns)
          if ($prev[0] eq 'sub') {
            $prev[1] .= $rseq[$pi]; $prev[4]++;
            $prev[2] .= $qseq[$pi]; $prev[6]++;
          } else {
            push(@blks, [@prev]);
            @prev = ('sub', $rseq[$pi], $qseq[$pi], $prev[4]+1, $prev[4]+1, $prev[6]+1, $prev[6]+1);
          }
        }
      }
      $blks[0][0] eq '' and shift(@blks);
      # Merge neighboring @blks:
      #   sub+ins+sub / ins+sub+ins => INS:ALN;
      #   sub+del+sub / del+sub+del => DEL:ALN;
      my @var_blks;
      for (my $j=0; $j<@blks; $j++) {
        $blks[$j][0] eq 'same' and next;
        if (scalar(@var_blks) == 0) {
          @var_blks = ([@{$blks[$j]}, $j, $blks[$j][0], $blks[$j][0]]); next;
        }
        if (!($j == $var_blks[-1][7]+1)) {
          push(@var_blks, [@{$blks[$j]}, $j, $blks[$j][0], $blks[$j][0]]); next;
        }
        if ($var_blks[-1][0] eq 'sub') {
          if ($blks[$j][0] eq 'ins') {
            $var_blks[-1][0] = 'ins';
            $var_blks[-1][2] .= $blks[$j][2];
            $var_blks[-1][6] = $blks[$j][6];
            $var_blks[-1][9] = $blks[$j][0];
          } elsif ($blks[$j][0] eq 'del') {
            $var_blks[-1][0] = 'del';
            $var_blks[-1][1] .= $blks[$j][1];
            $var_blks[-1][4] = $blks[$j][4];
            $var_blks[-1][9] = $blks[$j][0];
          } else {
            push(@var_blks, [@{$blks[$j]}, $j, $blks[$j][0], $blks[$j][0]]);
          }
        } elsif ($var_blks[-1][0] eq 'ins') {
          # This ins can be ins_only, ins_sub, sub_ins, ...
          if ($blks[$j][0] =~ m!^(sub|ins)$!i) {
            $var_blks[-1][1] .= $blks[$j][1];
            $var_blks[-1][2] .= $blks[$j][2];
            $var_blks[-1][4]  = $blks[$j][4];
            $var_blks[-1][6]  = $blks[$j][6];
            $var_blks[-1][9] = $blks[$j][0];
          } else {
            push(@var_blks, [@{$blks[$j]}, $j, $blks[$j][0], $blks[$j][0]]);
          }
        } elsif ($var_blks[-1][0] eq 'del') {
          # This del can be del_only, del_sub, sub_del, ...
          if ($blks[$j][0] =~ m!^(sub|del)$!i) {
            $var_blks[-1][1] .= $blks[$j][1];
            $var_blks[-1][2] .= $blks[$j][2];
            $var_blks[-1][4] = $blks[$j][4];
            $var_blks[-1][6] = $blks[$j][6];
            $var_blks[-1][9] = $blks[$j][0];
          } else {
            push(@var_blks, [@{$blks[$j]}, $j, $blks[$j][0], $blks[$j][0]]);
          }
        } else {
          push(@var_blks, [@{$blks[$j]}, $j, $blks[$j][0], $blks[$j][0]]);
        }
      }
      # Correct rS/rE/rSeq/qSeq/qS/qE.
      for my $v1 (@var_blks) {
        my $cor_rS = $rS + $v1->[3] - 1;
        my $cor_rE = $rS + $v1->[4] - 1;
        my $cor_qS = $qS + $v1->[5] - 1;
        my $cor_qE = $qS + $v1->[6] - 1;
        my $rev_qS = $qE - $v1->[6] + 1;
        my $rev_qE = $qE - $v1->[5] + 1;
        if      ($v1->[8] eq 'ins') {
          if ($v1->[1] eq '') {
            $v1->[1] = substr($refSeq{$refID}{'seq'}, $cor_rS-1, 1);
            $v1->[2] = $v1->[1] . $v1->[2];
            $v1->[3] = $cor_rS;   $v1->[4] = $cor_rS+1;
            $v1->[5] = $cor_qS;   $v1->[6] = $cor_qE;
            if ($strN == -1) {
              $v1->[5] = $rev_qS; $v1->[6] = $rev_qE;
            }
          } else {
            $v1->[3] = $cor_rS + 1; $v1->[4] = $cor_rE;
            $v1->[5] = $cor_qS;     $v1->[6] = $cor_qE;
            if ($strN == -1) {
              $v1->[5] = $rev_qS;   $v1->[6] = $rev_qE;
            }
          }
        } elsif ($v1->[8] eq 'del') {
          if ($v1->[2] eq '') {
            $v1->[2] = substr($refSeq{$refID}{'seq'}, $cor_rS-2, 1);
            $v1->[1] = $v1->[2] . $v1->[1];
            $v1->[3] = $cor_rS-1; $v1->[4] = $cor_rE;
            $v1->[5] = $cor_qS;   $v1->[6] = $cor_qS+1;
            if ($strN == -1) {
              $v1->[5] = $rev_qS-1; $v1->[6] = $rev_qE;
            }
          } else {
            $v1->[3] = $cor_rS;   $v1->[4] = $cor_rE;
            $v1->[5] = $cor_qS+1; $v1->[6] = $cor_qE;
            if ($strN == -1) {
             $v1->[5] = $rev_qS;  $v1->[6] = $rev_qE-1;
            }
          }
        } else {
          $v1->[3] = $cor_rS; $v1->[4] = $cor_rE;
          $v1->[5] = $cor_qS; $v1->[6] = $cor_qE;
          if ($strN == -1) {
            $v1->[5] = $rev_qS; $v1->[6] = $rev_qE;
          }
        }
      }
      # Output VCF record.
      for my $v1 (@var_blks) {
        my @v2=@$v1;
        my $svlen = length($v2[2]) - length($v2[1]);
        &boundary_inBad($refID, $v2[3], $v2[4], $qryID, $v2[5], $v2[6], $maxDist) and next;
        if ($v2[0] eq 'sub') {
          my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=UNK:ALN;END=$v2[4];SVLEN=$svlen;ALTPOS=$qryID:$v2[5]-$v2[6]";
          print STDOUT join("\t", $refID, $v2[3], '.', $v2[1], $v2[2], '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
        } elsif ($v2[0] eq 'ins') {
          my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=INS:ALN;END=$v2[4];SVLEN=$svlen;ALTPOS=$qryID:$v2[5]-$v2[6]";
          print STDOUT join("\t", $refID, $v2[3], '.', $v2[1], $v2[2], '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
        } elsif ($v1->[0] eq 'del') {
          my $infoTxt = "ALGORITHMS=NucDiff_struct;SVTYPE=DEL:ALN;END=$v2[4];SVLEN=$svlen;ALTPOS=$qryID:$v2[5]-$v2[6]";
          print STDOUT join("\t", $refID, $v2[3], '.', $v2[1], $v2[2], '.', 'PASS', $infoTxt, "GT", "1/1")."\n";
        } else {
          &stopErr("[Err] Undefined VAR type [$v1->[0]] in [@$v1]\n");
        }
      }
    }
  }
}

&fileSunhh::_rmtree($wdir);

#
# Subroutines.
#
sub filter_array {
  my ($ar, $hr) = @_;
  my @back;
  for (my $i=0; $i<@$ar; $i++) {
    defined $hr->{$i} and next;
    push(@back, $ar->[$i]);
  }
  return(@back);
}

sub sepTA8 {
  my %back;
  for (split(/;/, $_[0])) {
    m!^([^\s;]+)=([^\s;]+)$! or &stopErr("[Err] 1: $_\n");
    $back{$1} = $2;
  }
  defined $back{'ref_bases'}   and $back{'ref_bases'} = uc($back{'ref_bases'});
  defined $back{'query_bases'} and $back{'query_bases'} = uc($back{'query_bases'});
  return(\%back);
}# sepTA8

sub boundary_inBad {
  my ($rID, $rS, $rE, $qID, $qS, $qE, $dist) = @_;
  if (defined $bad_loc{'r'}{$rID}) {
    for my $a1 (@{$bad_loc{'r'}{$rID}}) {
      for my $b1 ($rS, $rE) {
        $a1->[0]-$dist <= $b1 and $a1->[0]+$dist >= $b1 and do {return 1;};
        $a1->[1]-$dist <= $b1 and $a1->[1]+$dist >= $b1 and do {return 1;};
      }
    }
  }
  if (defined $bad_loc{'q'}{$qID}) {
    for my $a1 (@{$bad_loc{'q'}{$qID}}) {
      for my $b1 ($qS, $qE) {
        $a1->[0]-$dist <= $b1 and $a1->[0]+$dist >= $b1 and do {return 1;};
        $a1->[1]-$dist <= $b1 and $a1->[1]+$dist >= $b1 and do {return 1;};
      }
    }
  }
  return(undef());
}# boundary_inBad

