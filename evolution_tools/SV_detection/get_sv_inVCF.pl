#!/usr/bin/perl
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
    if ($refLen == 1) {
      $altLen == 1 and do { $type //= 'SNP'; $alN2N{$j} = '.'; };
      if ($altLen-1 >= $minLen) {
        $type = 'long_SV'; $newAlN++; $alN2N{$j} = $newAlN;
      } else {
        (defined $type and $type =~ m!^(long_SV|long_MNP)$!) or $type = 'short_SV';
        $alN2N{$j} = '.';
      }
    } elsif ($altLen == 1) {
      if ($refLen-1 >= $minLen) {
        $type = 'long_SV'; $newAlN++; $alN2N{$j} = $newAlN;
      } else {
        (defined $type and $type =~ m!^(long_SV|long_MNP)$!) or $type = 'short_SV';
        $alN2N{$j} = '.';
      }
    } elsif ($refLen == $altLen) {
      if ($refLen >= $minLen) {
        (defined $type and $type eq 'long_SV') or $type = 'long_MNP';
        $newAlN++; $alN2N{$j} = $newAlN;
      } else {
        (defined $type and $type =~ m!^(long_SV|long_MNP)$!) or $type = 'short_MNP';
        $alN2N{$j} = '.';
      }
    } elsif ($refLen >= $minLen or $altLen >= $minLen) {
      $type = 'long_SV'; $newAlN++; $alN2N{$j} = $newAlN;
    } else {
      (defined $type and $type =~ m!^(long_SV|long_MNP)$!) or $type = 'short_SV';
      $alN2N{$j} = '.';
    }
  }# End for ()
  if ($type eq "long_SV" or $type eq "long_MNP") {
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

