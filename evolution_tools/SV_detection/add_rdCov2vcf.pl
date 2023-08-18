#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 min_read_depth_R min_read_depth_Q o_ref_struct.vcf Q2R.bam R2Q.bam > o_ref_struct.cov.vcf\n";

my ($minDR, $minDQ, $vcfFn, $rbamFn, $qbamFn) = @ARGV;
my $covtag = "RDCOVR${minDR}Q${minDQ}";

my $fh = &openFH($vcfFn);
my $has_header=0;
while (<$fh>) {
  if (m!^#!) {
    print STDOUT $_;
    m!^##INFO=\<ID=${covtag},! and $has_header=1;
    next;
  }
  if ($has_header == 0) {
    print STDOUT "##INFO=<ID=${covtag},Number=.,Type=Float,Description=\"Type of structural variant\">\n";
    $has_header = 1;
  }

  chomp;
  my @ta=split(/\t/, $_);
  # $ta[4] eq '<INV>' and do { print "$_\n"; next; };
  $ta[7] =~ m!END=(\d+)! or &stopErr("[Err] Bad input.vcf\n"); 
  my ($rID, $rS, $rE) = ($ta[0], $ta[1], $1);
  $ta[7] =~ m!ALTPOS=([^;\s]+):(\d+)\-(\d+)! or &stopErr("[Err] No ALTPOS: $ta[7]\n");
  my ($qID, $qS, $qE) = ($1, $2, $3);
  $ta[7] =~ m!SVTYPE=([^\s;]+)! or &stopErr("[Err] No SVTYPE: $ta[7]\n");
  my $svtype = $1;
  # my $refL = length($ta[3]);
  # my $altL = length($ta[4]);
  if ($svtype =~ m!^(DEL|DEL:ALN|DEL:REL|DEL:REL:UNK|DEL:REP|DEL:REP:UNK)$!i) {
    $rS ++;
    # In deletion, only $rbamFn matters.
    my $cR = &cnt_cov($rbamFn, $rID, $rS, $rE, $minDR);
    $ta[7] .= ";$covtag=$cR";
  } elsif ($svtype =~ m!^(INS|INS:ALN|INS:REL|INS:REL:UNK|INS:DUP|INS:DUP:UNK|INS:TAN)$!i) {
    # In insertion, only $qbamFn matters.
    my $cQ = &cnt_cov($qbamFn, $qID, $qS, $qE, $minDQ);
    $ta[7] .= ";$covtag=$cQ";
  } elsif ($svtype =~ m!^(UNK:SUB|DEL:SUB|INS:SUB|UNK:ALN|INV)$!i) {
    # In substitution and inversions, both $rbamFn and $qbamFn matter.
    my $cR = &cnt_cov($rbamFn, $rID, $rS, $rE, $minDR);
    my $cQ = &cnt_cov($qbamFn, $qID, $qS, $qE, $minDQ);
    $ta[7] .= ";$covtag=$cR,$cQ";
  } else {
    &stopErr("[Err] Unknown svtype [$svtype]\n");
  }
  print STDOUT join("\t", @ta)."\n";
}
close($fh);

sub cnt_cov {
  my ($bamFn, $id, $s, $e, $minD) = @_;
  open F1,'-|', "samtools depth $bamFn -r $id:$s-$e | cut -f 3" or &stopErr("[Err] CMD: $!\n");
  my $cov_len = 0;
  while (<F1>) {
    chomp;
    $_ >= $minD or next;
    $cov_len ++;
  }
  close F1;
  return(sprintf("%0.2f", 100*$cov_len/($e-$s+1)));
}
