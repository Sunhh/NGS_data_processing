#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use fileSunhh;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 ref.fa.gz qry.fa.gz sampleID 22CEXU12.aw_ref_snps.gff.gz > 22CEXU12.aw_ref_snps.vcf\n";

my $fnRFa = shift;
my $fnQFa = shift;
my $samID = shift;
my $fnGff = shift;

my %refSeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnRFa)};
for (keys %refSeq) { $refSeq{$_}{'seq'} =~ s!\s!!g; $refSeq{$_}{'len'}=length($refSeq{$_}{'seq'}); $refSeq{$_}{'seq'} = uc($refSeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded ref fa [$fnRFa]\n");

my %qrySeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnQFa)};
for (keys %qrySeq) { $qrySeq{$_}{'seq'} =~ s!\s!!g; $qrySeq{$_}{'len'}=length($qrySeq{$_}{'seq'}); $qrySeq{$_}{'seq'} = uc($qrySeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded qry fa [$fnQFa]\n");

my %errCnt;
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
HHHS

print STDOUT join("\t", "#CHROM", qw/POS ID REF ALT QUAL INFO FORMAT/, $samID)."\n";
while (<$fhGff>) {
  chomp;
  m!^\s*#! and next;
  my @inf = ("ALGORITHMS=NucDiff_snps");
  my @ta=split(/\t/, $_);
  my $h8 = &sepTA8($ta[8]);
  my (%refInf, %qryInf);
  $refInf{'chr'} = $ta[0];
  if      ($h8->{'Name'} =~ m!^deletion$!i) {
    push(@inf, "SVTYPE=DEL");
    $refInf{'pos'} = $ta[3]-1;
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, $ta[4]-$refInf{'pos'}+1);
    $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $h8->{'query_coord'}-1, 1);
    push(@inf, join("", "END=", $refInf{'pos'}+length($refInf{'base'})-1));
    if ($h8->{'query_dir'} == -1) {
      $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $h8->{'query_coord'}, 1);
      &fastaSunhh::rcSeq(\$qryInf{'base'}, 'rc');
    }
    unless ($qryInf{'base'} eq substr($refInf{'base'}, 0, 1)) {
      # I believe this doesn't happen a lot.
      # This happens when the AnchorWave makes wrong alignment. For example:
      # Aligned Ref: AAATGACTTAATTCAATAGAAT-----TGCTTTATATATATATATATATATATATATATTAAGTTACAATTAAAATA
      # Aligned Qry: AAATGACTTAATTCAATAGAATTGCTA----------------------------------AGTTACAATTAAAATA
      #   Resulting one deletion and one insertion.
      # True    Ref: AAATGACTTAATTCAATAGAATTGCT-----TTATATATATATATATATATATATATATTAAGTTACAATTAAAATA
      # True    Qry: AAATGACTTAATTCAATAGAATTGCT----------------------------------AAGTTACAATTAAAATA
      #   Resulting one deletion only.
      $errCnt{'err_del'} ++;
      if (length($refInf{'base'}) > 50) {
        $refInf{'base'} = substr($refInf{'base'}, 0, 49) . "...";
      }
      &tsmsg("[Wrn] problem line err_del:$errCnt{'err_del'}:ref [$refInf{'base'}], qry [$qryInf{'base'}]:\t$_\n");
      next;
    }
  } elsif ($h8->{'Name'} =~ m!^substitution$!i) {
    # There isn't a SVTYPE defined for equal-length substitution in VCF format.
    $refInf{'pos'} = $ta[3];
    $refInf{'base'} = $h8->{'ref_bases'};
    $qryInf{'base'} = $h8->{'query_bases'};
    ( defined $refInf{'base'} and $refInf{'base'} ne '' ) or &stopErr("[Err] No ref_bases found: $_\n");
    ( defined $qryInf{'base'} and $qryInf{'base'} ne '' ) or &stopErr("[Err] No query_bases found: $_\n");
    # query_dir=-1 doesn't change base information
  } elsif ($h8->{'Name'} =~ m!^insertion$!i) {
    push(@inf, "SVTYPE=INS");
    $refInf{'pos'} = $ta[3];
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, 1);
    $qryInf{'base'} = $refInf{'base'} . $h8->{query_bases};
    # RC alignment doesn't change the result.
  } elsif ($h8->{'Name'} =~ m!^(inserted_gap|gap)$!i) {
    # N-gap should be ignored.
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
  defined $back{'ref_bases'} and $back{'ref_bases'} = uc($back{'ref_bases'});
  defined $back{'query_bases'} and $back{'ref_bases'} = uc($back{'query_bases'});
  return(\%back);
}# sepTA8

