#!/usr/bin/perl
use strict;
use warnings;
use fileSunhh;
use LogInforSunhh;

!@ARGV and die "perl $0 readCovTag_RDCOVRnQn min_read_Cov_perc_60.0 cov_struct.vcf > filtered_struct.vcf\n";

my ($covtag, $minC, $vcfFn) = @ARGV;

my $fh = &openFH($vcfFn);
while (<$fh>) {
  if (m!^#!) {
    print STDOUT $_;
    next;
  }

  chomp;
  my @ta=split(/\t/, $_);
  $ta[4] eq '<INV>' and do { print STDOUT "$_\n"; next; };
  $ta[7] =~ m!${covtag}=([\d.,]+)! or do { print STDOUT "$_\n"; next; };
  my $txt1 = $1;
  if ($txt1 =~ m!^([\d.]+),([\d.]+)$!) {
    my ($rc, $qc) = ($1, $2);
    $rc > $minC and $qc > $minC and next;
    if ($rc > $minC) {
      length($ta[3]) == 1 and length($ta[4]) == 1 and next;
      $ta[3] = substr($ta[3], 0, 1);
      $ta[4] = $ta[3] . $ta[4];
      print STDOUT join("\t", @ta)."\n";
    } elsif ($qc > $minC) {
      length($ta[3]) == 1 and length($ta[4]) == 1 and next;
      $ta[4] = ".";
      print STDOUT join("\t", @ta)."\n";
    } else {
      print STDOUT "$_\n";
    }
  } elsif ($txt1 =~ m!^([\d.]+)$!) {
    my $cov = $1;
    if ($ta[7] =~ m!SVTYPE=(INS|DEL|UNK:ALN);! and length($ta[3]) < 10e3 and length($ta[4]) < 10e3 ) {
      print STDOUT "$_\n"; next;
    }
    $cov > $minC and next; # Some true short INS/DEL from tandem duplicates may be ignored!
    print STDOUT "$_\n";
  } else {
    &tsmsg("[Wrn] Error format of $covtag [$txt1]\n");
    print STDOUT "$_\n";
  }
}
close($fh);

