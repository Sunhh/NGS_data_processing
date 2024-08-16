#!/usr/bin/perl
# 4/19/2024: Update to include 1-to-2 ploidy change.
# Randomly choose two alleles from four alleles.
# CMD: bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP --threads 30 data/4x_to_tzn1hap.vcf.gz|perl tools/cnvt_ploidy_4to2.pl| bgzip -@ 30 -c > data/2x.vcf.gz
use strict;
use warnings;
use Getopt::Long;
use LogInforSunhh;

my %opts;
GetOptions(\%opts,
  "1to2!", # Ploidy 1to2.
  "4to2a!" # Default mode. Randomly choose two alleles.
);

my $htxt = <<HH;
################################################################################
bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP 4x.vcf.gz | perl $0 > 2x.vcf

 -1to2      [Boolean]
 -4to2a     [Boolean] default.
################################################################################
HH


-t and !@ARGV and &LogInforSunhh::usage($htxt);

$opts{'4to2a'} //= 1;

if ($opts{'1to2'}) {
  while (<>) {
    m!^#! and do {print; next;};
    chomp;
    my @ta=split(/\t/, $_);
    $ta[8] =~ m!^GT(:|$)! or die "Bad format: |$ta[8]|\n";
    $ta[8] = "GT";
    for my $b1 (@ta[9..$#ta]) {
      $b1 =~ m!^([\d.]+)(:|$)! or die "Bad genotype: |$b1| in line: $_\n";
      $b1 = "$1/$1";
    }
    print STDOUT join("\t", @ta)."\n";
  }
} elsif ($opts{'4to2a'}) {
  while (<>) {
    m!^#! and do {print; next;};
    chomp;
    my @ta=split(/\t/, $_);
    $ta[8] eq 'GT:AD:DP' or die "Bad FORMAT: |$ta[8]| $_\n";
    for my $b1 (@ta[9..$#ta]) {
      if ($b1 =~ s!^\./\./\./\.:!./.:! or $b1 =~ s!^0/0/0/0:!0/0:!) {
        next;
      }
      $b1 =~ s!^([\d\.]+)/([\d\.]+)/([\d\.]+)/([\d\.]+):!! or die "|$b1| in line:\n$_\n";
      my @c1=($1,$2,$3,$4);
      $b1 = $c1[int(rand(4))]."/".$c1[int(rand(4))] . ":$b1";
    }
    print STDOUT join("\t", @ta)."\n";
  }
}

