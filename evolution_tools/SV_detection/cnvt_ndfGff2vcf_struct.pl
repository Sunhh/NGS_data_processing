#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 ref.fa.gz qry.fa.gz sampleID 22CEXU12.aw_ref_struct.gff.gz > 22CEXU12.aw_ref_struct.vcf\n";

my $fnRFa = shift;
my $fnQFa = shift;
my $samID = shift;
my $fnGff = shift;

# Possible Names:
###  Accept: substitution;
###  Part  : gap/inserted_gap/unaligned_end/unaligned_beginning/...;
###  Accept: insertion/duplication/...;
###  Accept: deletion/collapsed_repeat/collapsed_tandem_repeat/...;
###  Ignore: relocation*/translocation*/inversion/reshuffling/...;

my %refSeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnRFa)};
for (keys %refSeq) { $refSeq{$_}{'seq'} =~ s!\s!!g; $refSeq{$_}{'len'}=length($refSeq{$_}{'seq'}); $refSeq{$_}{'seq'} = uc($refSeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded ref fa [$fnRFa]\n");

my %qrySeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnQFa)};
for (keys %qrySeq) { $qrySeq{$_}{'seq'} =~ s!\s!!g; $qrySeq{$_}{'len'}=length($qrySeq{$_}{'seq'}); $qrySeq{$_}{'seq'} = uc($qrySeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded qry fa [$fnQFa]\n");

my $fhGff = &openFH($fnGff);
print STDOUT <<'HHHS';
##fileformat=VCFv4.2
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=CNV,Description="Copy Number Polymorphism">
##ALT=<ID=CPX,Description="Complex SV">
##ALT=<ID=CTX,Description="Reciprocal chromosomal translocation">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">
##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">
##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">
##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">
##ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">
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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
HHHS

for (sort keys %refSeq) {
  print STDOUT "##contig=<ID=$_,length=$refSeq{$_}{'len'}>\n";
}

print STDOUT join("\t", "#CHROM", qw/POS ID REF ALT QUAL FILTER INFO FORMAT/, $samID)."\n";

while (<$fhGff>) {
  chomp;
  m!^\s*#! and next;
  my @ta=split(/\t/, $_);
  my @inf;
  push(@inf, "ALGORITHMS=NucDiff_struct");
  my $h8 = &sepTA8($ta[8]);
  my (%refInf, %qryInf);
  $refInf{'chr'} = $ta[0];
  if      ($h8->{'Name'} =~ m!^deletion$!i) {
    push(@inf, "SVTYPE=CPX;CPX_TYPE=CPX_TYPE_NDstruct_DEL");
    # [For ref_struct.gff] This is not a real deletion, but a residue from a big "substitution" that AnchorWave failed to align.
    $refInf{'pos'} = $ta[3]-1;
    $h8->{'query_coord'} =~ m!^\d+$! or &stopErr("[Err] Bad query_coord format for deletion [$h8->{'query_coord'}].\n");
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, $ta[4]-$refInf{'pos'}+1);
    $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $h8->{'query_coord'}-1, 1);
    if ($h8->{'query_dir'} == -1) {
      $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $h8->{'query_coord'}, 1);
      &fastaSunhh::rcSeq(\$qryInf{'base'}, 'rc');
    }
  } elsif ($h8->{'Name'} =~ m!^substitution$!i) {
    push(@inf, "SVTYPE=CPX;CPX_TYPE=CPX_TYPE_NDstruct_SUB");
    $refInf{'pos'} = $ta[3];
    $h8->{'query_coord'} =~ m!^(\d+)\-(\d+)$! or &stopErr("[Err] Bad query_coord format for substitution [$h8->{'query_coord'}].\n");
    my ($qs, $qe) = ($1, $2);
    $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $qs-1, $qe-$qs+1);
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, $ta[4]-$refInf{'pos'}+1);
    if ($h8->{'query_dir'} == -1) {
      &fastaSunhh::rcSeq(\$qryInf{'base'}, 'rc');
    }
  } elsif ($h8->{'Name'} =~ m!^insertion$!i) {
    push(@inf, "SVTYPE=CPX;CPX_TYPE=CPX_TYPE_NDstruct_INS");
    $refInf{'pos'} = $ta[3];
    $h8->{'query_coord'} =~ m!^(\d+)\-(\d+)$! or &stopErr("[Err] Bad query_coord format for insertion [$h8->{'query_coord'}].\n");
    my ($qs, $qe) = ($1, $2);
    $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $qs-1-1, $qe-$qs+1+1);
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, 1);
    if ($h8->{'query_dir'} == -1) {
      $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $qs-1, $qe-$qs+1+1);
      &fastaSunhh::rcSeq(\$qryInf{'base'}, 'rc');
    }
  } elsif ($h8->{'Name'} =~ m!^(unaligned_beginning|unaligned_end)$!i) {
    $h8->{'query_dir'} == 1 or &stopErr("[Err] query_dir is not 1 at $h8->{'Name'}\n");
    push(@inf, "SVTYPE=CPX;CPX_TYPE=CPX_TYPE_NDstruct_uALN");
    $refInf{'pos'} = $ta[3];
    $h8->{'query_coord'} =~ m!^(\d+)\-(\d+)$! or &stopErr("[Err] Bad query_coord format for unaligned_beginning/unaligned_end [$h8->{'query_coord'}].\n");
    my ($qs, $qe) = ($1, $2);
    $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $qs-1, $qe-$qs+1);
    if      ($h8->{'Name'} =~ m!^(unaligned_beginning)$!i) {
      $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, 0, $ta[4]);
    } elsif ($h8->{'Name'} =~ m!^(unaligned_end)$!i) {
      $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $ta[3]-1);
    } else {
      &stopErr("[Err] Why here 1!\n");
    }
  } elsif ($h8->{'Name'} =~ m!^(inversion|relocation.*|translocation.*|reshuffling.*|gap|inserted_gap|duplication|collapsed_repeat|collapsed_tandem_repeat)$!i) {
    # N-gap should be ignored.
    # duplication: https://github.com/uio-cels/NucDiff/wiki/Duplications#duplications
    #   There aren't many and are too complicated to introduce new variant sequence, so I ignore them.
    # collapsed_repeat: https://github.com/uio-cels/NucDiff/wiki/Collapsed-repeats
    #   Similar reason as duplication.
    # collapsed_tandem_repeat:
    #   case 1: https://github.com/uio-cels/NucDiff/wiki/Collapsed-tandem-repeats-(case-1)
    #   case 2: https://github.com/uio-cels/NucDiff/wiki/Collapsed-tandem-repeats-(case-2)
    #   Too complicated.
    next;
  } else {
    &stopErr("[Err] Wrong name [$h8->{'Name'}]\n");
  }
  for my $a1 (qw/chr pos base/) {
    defined $refInf{$a1} or die "aaa: $a1: $_\n";
  }
  defined $qryInf{'base'} or die "q_base: $_\n";
  scalar(@inf) == 0 and @inf = (".");
  print STDOUT join("\t", $refInf{'chr'}, $refInf{'pos'}, '.', $refInf{'base'}, $qryInf{'base'}, '.', '.', join(";", @inf), 'GT', '1/1')."\n";
}
close ($fhGff);


sub sepTA8 {
  my %back;
  for (split(/;/, $_[0])) {
    m!^([^\s;]+)=([^\s;]+)$! or &stopErr("[Err] 1: $_\n");
    $back{$1} = $2;
  }
  defined $back{'ref_bases'}   and $back{'ref_bases'} = uc($back{'ref_bases'});
  defined $back{'query_bases'} and $back{'ref_bases'} = uc($back{'query_bases'});
  return(\%back);
}# sepTA8

