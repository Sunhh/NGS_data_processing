#!/usr/bin/perl
# [3/23/2022] From liftoff output gff3 file, retrieve old mRNA ID, new mRNA ID, coverage_by_liftoff, cds length on target genome, strand on target genome.
# [10/11/2022] Allow other sections between "ID" and "Parent" in gff3 $ta[8] field.
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

-t and !@ARGV and die "perl $0 output/mapCDS.ASM1003_PI537277.to.ASM1002_97103.gff3 > output/mapCDS.ASM1003_PI537277.to.ASM1002_97103.gff3.tbl\n";

my %skip_feat;

my (%g, %m, %c);
while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  if ($ta[2] =~ m!^gene$!i) {
    $ta[8] =~ m!^ID=([^\s;]+);.*coverage=([\d.]+);sequence_ID=([\d.]+);.*extra_copy_number=(\d+)(?:;|$)! or die "bad gene ta8: $_\n";
    my ($tgt_geneID, $lo_cov, $lo_ident, $lo_cpN) = ($1, $2, $3, $4);
    my $ref_geneID = $tgt_geneID;
    if ($lo_cpN > 0) {
      $ref_geneID =~ s!_$lo_cpN$!! or die "Failed to parse 1: $ta[8]\n";
    }
    defined $g{$tgt_geneID} and die "repeat G $tgt_geneID: $ta[8]\n";
    $g{$tgt_geneID} = [$ref_geneID, $lo_cov, $lo_ident, $lo_cpN];
  } elsif ($ta[2] =~ m!^mRNA$!i) {
    $ta[8] =~ m!^ID=([^\s;]+);(?:[^\s;]+=[^\s;]+;)*Parent=([^\s;]+);.*extra_copy_number=(\d+)(?:;|$)! or die "bad mrna ta8: $_\n";
    my ($tgt_mrnaID, $tgt_geneID, $lo_cpN) = ($1, $2, $3);
    my $ref_mrnaID = $tgt_mrnaID;
    if ($lo_cpN > 0) {
      $ref_mrnaID =~ s!_$lo_cpN$!! or die "Failed to parse 2: $ta[8]\n";
    }
    defined $m{$tgt_mrnaID} and die "repeat M $tgt_mrnaID: $ta[8]\n";
    $m{$tgt_mrnaID} = [$ref_mrnaID, $tgt_geneID, $lo_cpN, @ta[0,3,4,6]];
  } elsif ($ta[2] =~ m!^CDS$!i) {
    $ta[8] =~ m!(?:;|^)\s*Parent=([^\s;]+)! or die "bad cds ta8: $_\n";
    my $tgt_mrnaID = $1;
    $c{$tgt_mrnaID}[0] += $ta[4]-$ta[3]+1;
    $c{$tgt_mrnaID}[1] ++;
    push(@{$c{$tgt_mrnaID}[2]}, [@ta[3,4]]);
  } else {
    if (!defined $skip_feat{$ta[2]}) {
      &tsmsg("[Wrn] Skip feature [$ta[2]]\n");
      $skip_feat{$ta[2]} = 1;
    }
  }
}
print STDOUT join("\t", qw/new_mrnaID coverage identity new_cdsLen new_exonNum ref_mrnaID ref_geneID exCopyN new_chrID new_chrS new_chrE new_chrStr new_CDSBlocks/)."\n";
for my $tgt_mID (sort { $m{$a}[3] cmp $m{$b}[3] || $m{$a}[4] <=> $m{$b}[4] || $m{$a}[5] <=> $m{$b}[5] } keys %m) {
  my ($ref_mID, $tgt_gID, $lo_cpN, $tgt_chrID, $tgt_chrS, $tgt_chrE, $tgt_chrStr) = @{$m{$tgt_mID}};
  defined $g{$tgt_gID} or die "Failed to find gene for $tgt_mID\n";
  my ($ref_gID, $lo_cov, $lo_ident, $lo_cpN_g) = @{$g{$tgt_gID}};
  $lo_cpN == $lo_cpN_g or die "something wrong 1\n";
  my $txt_cdsBlks;
  if ($tgt_chrStr eq '-') {
    $txt_cdsBlks = join(";", map { "$_->[0],$_->[1]" } sort { $b->[1] <=> $a->[1] || $b->[0] <=> $a->[0] } @{$c{$tgt_mID}[2]});
  } else {
    $txt_cdsBlks = join(";", map { "$_->[0],$_->[1]" } sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$c{$tgt_mID}[2]});
  }
  print STDOUT join("\t", $tgt_mID, $lo_cov, $lo_ident, $c{$tgt_mID}[0], $c{$tgt_mID}[1], $ref_mID, $ref_gID, $lo_cpN, $tgt_chrID, $tgt_chrS, $tgt_chrE, $tgt_chrStr, $txt_cdsBlks)."\n";
}


