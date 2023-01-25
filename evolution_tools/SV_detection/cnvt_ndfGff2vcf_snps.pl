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
for (keys %refSeq) { $refSeq{$_}{'seq'} =~ s!\s!!g; $refSeq{$_}{'len'}=length($refSeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded ref fa [$fnRFa]\n");

my %qrySeq = %{$fs_obj->save_seq_to_hash('faFile' => $fnQFa)};
for (keys %qrySeq) { $qrySeq{$_}{'seq'} =~ s!\s!!g; $qrySeq{$_}{'len'}=length($qrySeq{$_}{'seq'}); }
# &tsmsg("[Msg] Loaded qry fa [$fnQFa]\n");

my %errCnt;
my $fhGff = &openFH($fnGff);
print STDOUT join("\t", "#CHROM", qw/POS ID REF ALT QUAL INFO FORMAT/, $samID)."\n";
while (<$fhGff>) {
  chomp;
  m!^\s*#! and next;
  my @ta=split(/\t/, $_);
  my $h8 = &sepTA8($ta[8]);
  my (%refInf, %qryInf);
  $refInf{'chr'} = $ta[0];
  if      ($h8->{'Name'} =~ m!^deletion$!i) {
    $refInf{'pos'} = $ta[3]-1;
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, $ta[4]-$refInf{'pos'}+1);
    $qryInf{'base'} = substr($qrySeq{ $h8->{'query_sequence'} }{'seq'}, $h8->{'query_coord'}-1, 1);
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
    $refInf{'pos'} = $ta[3];
    $refInf{'base'} = $h8->{'ref_bases'};
    $qryInf{'base'} = $h8->{'query_bases'};
    # query_dir=-1 doesn't change base information
  } elsif ($h8->{'Name'} =~ m!^insertion$!i) {
    $refInf{'pos'} = $ta[3];
    $refInf{'base'} = substr($refSeq{ $refInf{'chr'} }{'seq'}, $refInf{'pos'}-1, 1);
    $qryInf{'base'} = $refInf{'base'} . $h8->{query_bases};
    # RC alignment doesn't change the result.
  } elsif ($h8->{'Name'} =~ m!^(inserted_gap|gap)$!i) {
    # N-gap should be ignored.
    next;
  } else {
    # edit here.
    &stopErr("[Err] Wrong name [$h8->{'Name'}]\n");
  }
  for my $a1 (qw/chr pos base/) {
    defined $refInf{$a1} or die "aaa: $a1: $_\n";
  }
  defined $qryInf{'base'} or die "q_base: $_\n";
  print STDOUT join("\t", $refInf{'chr'}, $refInf{'pos'}, '.', $refInf{'base'}, $qryInf{'base'}, '.', '.', '.', 'GT', '1/1')."\n";
}
close ($fhGff);


sub sepTA8 {
  my %back;
  for (split(/;/, $_[0])) {
    m!^([^\s;]+)=([^\s;]+)$! or &stopErr("[Err] 1: $_\n");
    $back{$1} = $2;
  }
  return(\%back);
}# sepTA8

