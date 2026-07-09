#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $FLOOR = 0.01;

sub usage {
	my $msg = shift // '';
	print STDERR $msg if $msg;
	print STDERR <<"HELP";

Usage:
  perl deg_call.pl --tpm <TPM> --pair <comp_grpPair_list> --fdr <fdr.list> [--out <file>] [--add_repTPM]

Reads per-replicate TPM values, a list of DESeq2 comparison pairs, and a DESeq2
FDR table, then writes a DEG-call table (one row per gene).

Required parameters:
  --tpm   <file>   TPM table. Row1 = header; col1 = GeneID; other columns are
                   per-replicate TPMs named <sample>_rep<N> (e.g. LA1357_D1_C_rep1).
  --pair  <file>   Comparison list. Tab-separated; col2 = sample1, col3 = sample2.
                   A sample name is the replicate name minus the trailing _rep<N>.
  --fdr   <file>   FDR table. Row1 = header; col1 = GeneID; comparison columns are
                   named ds.<sample1>_VS_<sample2>. Values may be NA.

Optional parameters:
  --out   <file>   Write output here (default: STDOUT).
  --add_repTPM     For each comparison, append every individual replicate TPM
                   (sample1 reps then sample2 reps) after that comparison's DEG
                   column, so replicate TPMs follow their DEG/FDR results.
  --help           Show this message and exit.

Output columns (per comparison, sample1 vs sample2):
  <cmp>.v1    mean TPM of sample1 replicates
  <cmp>.v2    mean TPM of sample2 replicates
  <cmp>.R     v2/v1, with v1 and v2 each floored to a minimum of $FLOOR first
  <cmp>.FDR   FDR from column ds.<sample1>_VS_<sample2>
  <cmp>.DEG   up   if R > 2   & FDR < 0.05
              down if R < 0.5 & FDR < 0.05
              not  otherwise (FDR is NA, or FDR >= 0.05, or 0.5 <= R <= 2)
  (with --add_repTPM: followed by each replicate TPM, column-named as in --tpm)

HELP
	exit($msg ? 1 : 0);
}

my ($tpm_file, $comp_file, $fdr_file, $out_file, $add_repTPM, $help);
GetOptions(
	'tpm=s'       => \$tpm_file,
	'pair=s'      => \$comp_file,
	'fdr=s'       => \$fdr_file,
	'out=s'       => \$out_file,
	'add_repTPM!' => \$add_repTPM,
	'help|h'      => \$help,
) or usage("ERROR: bad command line.\n");

usage() if $help;
usage("ERROR: --tpm is required.\n")  unless defined $tpm_file;
usage("ERROR: --pair is required.\n") unless defined $comp_file;
usage("ERROR: --fdr is required.\n")  unless defined $fdr_file;

# ---- 1. Read TPM header: map sample name to replicate column indices ----
open my $TPM, '<', $tpm_file or die "Cannot open --tpm $tpm_file: $!\n";
my $tpm_head = <$TPM>;
chomp $tpm_head;
my @tpm_cols = split /\t/, $tpm_head;
my %sample2idx;                     # sample name => [col indices, in file order]
for my $i (1 .. $#tpm_cols) {
	(my $sample = $tpm_cols[$i]) =~ s/_rep\d+$//;
	push @{ $sample2idx{$sample} }, $i;
}

# gene => sample => mean TPM
my %tpm_mean;
while (my $line = <$TPM>) {
	chomp $line;
	my @f = split /\t/, $line;
	my $gid = $f[0];
	for my $sample (keys %sample2idx) {
		my $sum = 0;
		my $n   = 0;
		for my $i (@{ $sample2idx{$sample} }) {
			my $v = $f[$i];
			next unless defined $v && $v ne '' && $v ne 'NA';
			$sum += $v;
			$n++;
		}
		$tpm_mean{$gid}{$sample} = $n ? $sum / $n : 'NA';
	}
}
close $TPM;

# ---- 2. Read FDR list header + values ----
open my $FDR, '<', $fdr_file or die "Cannot open --fdr $fdr_file: $!\n";
my $fdr_head = <$FDR>;
chomp $fdr_head;
my @fdr_cols = split /\t/, $fdr_head;
my %fdr_colidx;                     # column name => index
for my $i (1 .. $#fdr_cols) {
	$fdr_colidx{ $fdr_cols[$i] } = $i;
}

my %fdr_val;                        # gene => [fields]
while (my $line = <$FDR>) {
	chomp $line;
	my @f = split /\t/, $line;
	my $gid = $f[0];
	$fdr_val{$gid} = \@f;
}
close $FDR;

# ---- 3. Read comparison list ----
open my $CMP, '<', $comp_file or die "Cannot open --pair $comp_file: $!\n";
my @comps;                          # each: [sample1, sample2, fdr_colname]
while (my $line = <$CMP>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @f = split /\t/, $line;
	my ($s1, $s2) = ($f[1], $f[2]);
	my $fdr_name = "ds.${s1}_VS_${s2}";
	exists $fdr_colidx{$fdr_name}
		or warn "WARN: FDR column '$fdr_name' not found in --fdr $fdr_file\n";
	exists $sample2idx{$s1}
		or warn "WARN: sample '$s1' has no replicates in --tpm $tpm_file\n";
	exists $sample2idx{$s2}
		or warn "WARN: sample '$s2' has no replicates in --tpm $tpm_file\n";
	push @comps, [$s1, $s2, $fdr_name];
}
close $CMP;

# ---- 4. Open output ----
my $OUT;
if (defined $out_file) {
	open $OUT, '>', $out_file or die "Cannot open --out $out_file: $!\n";
} else {
	$OUT = \*STDOUT;
}

# ---- 5. Header ----
my @head = ('GeneID');
for my $c (@comps) {
	my ($s1, $s2) = ($c->[0], $c->[1]);
	my $tag = "${s1}_VS_${s2}";
	push @head, "${tag}.v1", "${tag}.v2", "${tag}.R", "${tag}.FDR", "${tag}.DEG";
	if ($add_repTPM) {
		push @head, map { $tpm_cols[$_] } @{ $sample2idx{$s1} } if exists $sample2idx{$s1};
		push @head, map { $tpm_cols[$_] } @{ $sample2idx{$s2} } if exists $sample2idx{$s2};
	}
}
print {$OUT} join("\t", @head), "\n";

# ---- 6. One row per gene (TPM file order) ----
open $TPM, '<', $tpm_file or die "Cannot open --tpm $tpm_file: $!\n";
<$TPM>;                             # skip header
while (my $line = <$TPM>) {
	chomp $line;
	my @f = split /\t/, $line;
	my $gid = $f[0];
	my @out = ($gid);
	for my $c (@comps) {
		my ($s1, $s2, $fdr_name) = @$c;

		my $v1 = $tpm_mean{$gid}{$s1};
		my $v2 = $tpm_mean{$gid}{$s2};
		$v1 = 'NA' unless defined $v1;
		$v2 = 'NA' unless defined $v2;

		# floored values for the ratio
		my $R = 'NA';
		if ($v1 ne 'NA' && $v2 ne 'NA') {
			my $v1f = $v1 < $FLOOR ? $FLOOR : $v1;
			my $v2f = $v2 < $FLOOR ? $FLOOR : $v2;
			$R = $v2f / $v1f;
		}

		# FDR lookup
		my $fdr = 'NA';
		if (exists $fdr_colidx{$fdr_name} && defined $fdr_val{$gid}) {
			my $idx = $fdr_colidx{$fdr_name};
			my $raw = $fdr_val{$gid}[$idx];
			$fdr = (defined $raw && $raw ne '') ? $raw : 'NA';
		}

		# DEG call
		my $deg = 'not';
		if ($fdr ne 'NA' && $fdr < 0.05 && $R ne 'NA') {
			if    ($R > 2)   { $deg = 'up';   }
			elsif ($R < 0.5) { $deg = 'down'; }
		}

		push @out, $v1, $v2, $R, $fdr, $deg;

		if ($add_repTPM) {
			push @out, map { defined $f[$_] ? $f[$_] : 'NA' } @{ $sample2idx{$s1} } if exists $sample2idx{$s1};
			push @out, map { defined $f[$_] ? $f[$_] : 'NA' } @{ $sample2idx{$s2} } if exists $sample2idx{$s2};
		}
	}
	print {$OUT} join("\t", @out), "\n";
}
close $TPM;
close $OUT if defined $out_file;
