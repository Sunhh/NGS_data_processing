#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use Getopt::Long;

# Keep only polymorphic sites (>= 2 distinct non-N genotypes among the genotype columns).
# Merges the old rm_same_site.pl + rm_same_site_hete2N.pl (they differed only in whether
# heterozygous codes are collapsed to N first) into one, selectable with -hete2N.

my %opts;
GetOptions(\%opts, "hete2N!", "geno_col:i", "help!");
my $geno_col = defined $opts{'geno_col'} ? $opts{'geno_col'} : 2;

my $help = <<HH;

perl $0 [-hete2N] [-geno_col $geno_col] in_pasted_1col.tbl > snp.tbl

  Keeps only sites with >= 2 distinct non-N genotypes from column geno_col onward.
  -hete2N     Collapse any genotype that is not a single ATGC base or an indel (*, +) to
              'N' first (i.e. heterozygous codes become N; indels accepted).
              Default: compare genotypes as-is.
  -geno_col   First genotype column, 0-based [$geno_col].

HH
($opts{'help'} or (-t STDIN and !@ARGV)) and &LogInforSunhh::usage($help);

while (<>) {
	$. % 1e6 == 1 and &tsmsg("[Msg] Reading $. lines.\n");
	s/[^\S\t]+$//;
	my @ta = split(/\t/, $_);
	if ( $ta[0] eq 'chr' ) {
		print "$_\n";
		next;
	}
	my $base = 'N';
	my $has_diff = 0;
	for (my $i=$geno_col; $i<@ta; $i++) {
		$ta[$i] = uc($ta[$i]);
		$opts{'hete2N'} and ( $ta[$i] =~ m/^[ATGC]$|\*|\+/ or $ta[$i] = 'N' );
		$ta[$i] eq 'N' and next;
		$base eq 'N' and $base = $ta[$i];
		$base ne $ta[$i] and do { $has_diff = 1; last; };
	}
	$has_diff == 1 and print "$_\n";
}
