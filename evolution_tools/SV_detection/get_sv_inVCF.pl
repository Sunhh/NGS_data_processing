#!/usr/bin/perl
# 3/10/2023: Better count SV size.
use strict;
use warnings;

!@ARGV and die "perl $0 20_minLength allele_merged.vcf > good_SVs.vcf\n";

my $minLen = shift;

while (<>) {
  if (m!^\s*(#|$)!) {
    print STDOUT $_;
    next;
  }
  chomp;
  my @ta=split(/\t/, $_);
  $ta[3] = uc($ta[3]);
  $ta[4] = uc($ta[4]);
  my $refA = $ta[3];
  my $refLen = length($ta[3]);
  my $type;
  my @altAlleles = split(/,/, $ta[4]);
  my %alN2N; $alN2N{'0'} = 0;
  my $newAlN = 0;
  for (my $i=0; $i<@altAlleles; $i++) {
    my $j = $i+1;
    my $altA = $altAlleles[$i];
    my $altLen = length($altA);
    if (index($refA, $altA)==0) {
      if ($refLen-$altLen >= $minLen) {
        $type = 'good'; $newAlN++; $alN2N{$j} = $newAlN;
      } else {
        (defined $type and $type =~ m!^good$!) or $type = 'bad'; $alN2N{$j} = '.';
      }
    } elsif (index($altA, $refA)==0) {
      if ($altLen-$refLen >= $minLen) {
        $type = 'good'; $newAlN++; $alN2N{$j} = $newAlN;
      } else {
        (defined $type and $type =~ m!^good$!) or $type = 'bad'; $alN2N{$j} = '.';
      }
    } elsif ($refLen >= $minLen or $altLen >= $minLen) {
      $type = 'good'; $newAlN++; $alN2N{$j} = $newAlN; # long_MNP is also OK here.
    } else {
      # 'short_MNP' and 'SNP' can here too.
      (defined $type and $type =~ m!^good$!) or $type = 'bad'; $alN2N{$j} = '.';
    }
  }# End for ()
  if ($type eq "good") {
    if ($newAlN < scalar(@altAlleles)) {
      my @tb;
      for (my $i=0; $i<@altAlleles; $i++) {
        my $j=$i+1;
        $alN2N{"$j"} eq "." and next;
        push(@tb, $altAlleles[$i]);
      }
      $ta[4] = join(",", @tb);
      for (my $i=9; $i<@ta; $i++) {
        $ta[$i] eq './.' and next;
        $ta[$i] =~ m!^(\d+)\/(\d+)$! or die "[Err] Bad allele [$ta[$i]] in line: $_\n";
        if ($alN2N{$1} eq "." or $alN2N{$2} eq ".") {
          $ta[$i] = './.';
        } else {
          $ta[$i] = $alN2N{$1}. "/". $alN2N{$2};
        }
      }
    }
    print STDOUT join("\t", @ta)."\n";
  }
}

