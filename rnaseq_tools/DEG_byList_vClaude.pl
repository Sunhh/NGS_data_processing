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
                   per-replicate TPMs named <group>_rep<N> / <group>_Rep<N> (case-insensitive; e.g. S1_Rep1).
  --pair  <file>   Comparison list. Tab-separated; col2 = group1, col3 = group2.
                   A group is the replicate column name minus the trailing _rep<N> (case-insensitive).
  --fdr   <file>   FDR table. Row1 = header; col1 = GeneID; each comparison is one
                   column named <method>.<group1>_VS_<group2> (ds. DESeq2, et. edgeR
                   classic, eg. edgeR GLM; a bare <group1>_VS_<group2> also works).
                   Values may be NA.

Optional parameters:
  --out   <file>   Write output here (default: STDOUT).
  --add_repTPM     For each comparison, append every individual replicate TPM
                   (group1 reps then group2 reps) after that comparison's DEG
                   column, so replicate TPMs follow their DEG/FDR results.
  --fdr_cut <f>    FDR significance cutoff for a DEG call (default: 0.05).
  --fc      <f>    Minimum TPM-mean fold change for a DEG call (default: 2);
                   up needs R > fc, down needs R < 1/fc.
  --floor   <f>    Floor applied to each mean TPM before the ratio (default: $FLOOR).
  --help           Show this message and exit.

Output columns (per comparison, group1 vs group2):
  <cmp>.v1    mean TPM of group1 replicates
  <cmp>.v2    mean TPM of group2 replicates
  <cmp>.R     v2/v1, with v1 and v2 each floored to a minimum of $FLOOR first
  <cmp>.FDR   FDR from column ds.<group1>_VS_<group2>
  <cmp>.DEG   up   if R > fc     & FDR < fdr_cut
              down if R < 1/fc   & FDR < fdr_cut
              not  otherwise (FDR is NA, or FDR >= fdr_cut, or 1/fc <= R <= fc)
  (with --add_repTPM: followed by each replicate TPM, column-named as in --tpm)

HELP
	exit($msg ? 1 : 0);
}

my ($tpm_file, $comp_file, $fdr_file, $out_file, $add_repTPM, $help);
my ($fdr_cut, $fc, $floor);
GetOptions(
	'tpm=s'       => \$tpm_file,
	'pair=s'      => \$comp_file,
	'fdr=s'       => \$fdr_file,
	'out=s'       => \$out_file,
	'add_repTPM!' => \$add_repTPM,
	'fdr_cut=f'   => \$fdr_cut,
	'fc=f'        => \$fc,
	'floor=f'     => \$floor,
	'help|h'      => \$help,
) or usage("ERROR: bad command line.\n");
$fdr_cut //= 0.05;
$fc      //= 2;
$floor   //= $FLOOR;

usage() if $help;
usage("ERROR: --tpm is required.\n")  unless defined $tpm_file;
usage("ERROR: --pair is required.\n") unless defined $comp_file;
usage("ERROR: --fdr is required.\n")  unless defined $fdr_file;

# ---- 1. Read TPM header: map group name (replicate name minus _rep<N>) to replicate column indices ----
open my $TPM, '<', $tpm_file or die "Cannot open --tpm $tpm_file: $!\n";
my $tpm_head = <$TPM>;
chomp $tpm_head;
my @tpm_cols = split /\t/, $tpm_head;
my %group2idx;                     # group name (replicate name minus _rep<N>) => [col indices, in file order]
for my $i (1 .. $#tpm_cols) {
	(my $group = $tpm_cols[$i]) =~ s/_rep\d+$//i;   # case-insensitive: _rep1 or _Rep1
	push @{ $group2idx{$group} }, $i;
}

# gene => sample => mean TPM
my %tpm_mean;
while (my $line = <$TPM>) {
	chomp $line;
	my @f = split /\t/, $line;
	my $gid = $f[0];
	for my $group (keys %group2idx) {
		my $sum = 0;
		my $n   = 0;
		for my $i (@{ $group2idx{$group} }) {
			my $v = $f[$i];
			next unless defined $v && $v ne '' && $v ne 'NA';
			$sum += $v;
			$n++;
		}
		$tpm_mean{$gid}{$group} = $n ? $sum / $n : 'NA';
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
my @comps;                          # each: [group1, group2, fdr_colname]
while (my $line = <$CMP>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my @f = split /\t/, $line;
	my ($g1, $g2) = ($f[1], $f[2]);
	my $suffix = "${g1}_VS_${g2}";
	# Match this pair's FDR column, tolerating any method prefix (ds./et./eg./none).
	my $fdr_name;
	for my $cand ("ds.$suffix", "et.$suffix", "eg.$suffix", $suffix) {
		if (exists $fdr_colidx{$cand}) { $fdr_name = $cand; last; }
	}
	unless (defined $fdr_name) {
		my @m = grep { $_ eq $suffix || /\.\Q$suffix\E$/ } keys %fdr_colidx;
		if (@m) { $fdr_name = $m[0]; warn "WARN: multiple FDR columns match '$suffix' (@m); using $fdr_name\n" if @m > 1; }
	}
	unless (defined $fdr_name) {
		$fdr_name = "ds.$suffix";
		warn "WARN: FDR column for '$suffix' not found in --fdr $fdr_file\n";
	}
	exists $group2idx{$g1}
		or warn "WARN: group '$g1' has no replicates in --tpm $tpm_file\n";
	exists $group2idx{$g2}
		or warn "WARN: group '$g2' has no replicates in --tpm $tpm_file\n";
	push @comps, [$g1, $g2, $fdr_name];
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
	my ($g1, $g2) = ($c->[0], $c->[1]);
	my $tag = "${g1}_VS_${g2}";
	push @head, "${tag}.v1", "${tag}.v2", "${tag}.R", "${tag}.FDR", "${tag}.DEG";
	if ($add_repTPM) {
		push @head, map { $tpm_cols[$_] } @{ $group2idx{$g1} } if exists $group2idx{$g1};
		push @head, map { $tpm_cols[$_] } @{ $group2idx{$g2} } if exists $group2idx{$g2};
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
		my ($g1, $g2, $fdr_name) = @$c;

		my $v1 = $tpm_mean{$gid}{$g1};
		my $v2 = $tpm_mean{$gid}{$g2};
		$v1 = 'NA' unless defined $v1;
		$v2 = 'NA' unless defined $v2;

		# floored values for the ratio
		my $R = 'NA';
		if ($v1 ne 'NA' && $v2 ne 'NA') {
			my $v1f = $v1 < $floor ? $floor : $v1;
			my $v2f = $v2 < $floor ? $floor : $v2;
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
		if ($fdr ne 'NA' && $fdr < $fdr_cut && $R ne 'NA') {
			if    ($R > $fc)     { $deg = 'up';   }
			elsif ($R < 1/$fc)   { $deg = 'down'; }
		}

		push @out, $v1, $v2, $R, $fdr, $deg;

		if ($add_repTPM) {
			push @out, map { defined $f[$_] ? $f[$_] : 'NA' } @{ $group2idx{$g1} } if exists $group2idx{$g1};
			push @out, map { defined $f[$_] ? $f[$_] : 'NA' } @{ $group2idx{$g2} } if exists $group2idx{$g2};
		}
	}
	print {$OUT} join("\t", @out), "\n";
}
close $TPM;
close $OUT if defined $out_file;
