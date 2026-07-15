#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Unified read-count -> normalized-expression converter (TPM and/or RPKM).
# Supersedes/merges cnvt_cnt_to_TPM.pl (Perl, long-format inputs) and
# cnvt_featureCounts_to_tpm.r (R, featureCounts-only, TPM-only).

my %o;
GetOptions(\%o,
  'cntFn=s',      # input count table (required)
  'fmt=s',        # auto | featurecounts | matrix   (default auto)
  'lenFn=s',      # gene-length file (required for --fmt matrix)
  'libFn=s',      # optional: per-sample true library size (total mapped reads)
  'out_type=s',   # tpm | rpkm | both               (default tpm)
  'outFn=s',      # primary output (default STDOUT)
  'rpkm_out=s',   # RPKM output when --out_type both
  'skipEleID=s',  # optional file: extra gene IDs to drop
  'help!',
);

sub usage {
  print STDERR <<"HELP";

Convert a gene read-count matrix to normalized expression (TPM and/or RPKM).

Usage:
  perl $0 --cntFn <counts> [--fmt auto] [--lenFn len] [--libFn lib] \\
          [--out_type tpm|rpkm|both] [--outFn out] [--rpkm_out rpkm.tsv]

Input count table (--cntFn), auto-detected (or force with --fmt):
  featurecounts  featureCounts output: columns Geneid,Chr,Start,End,Strand,Length,<sample..>.
                 Gene length is read from the Length column; sample counts are columns 7+.
  matrix         plain matrix: col1 = GeneID, other columns = per-sample counts.
                 Requires --lenFn for gene lengths.

Gene length (--lenFn; needed for --fmt matrix, ignored for featurecounts):
  2 columns:  GeneID <TAB> length
  or 3-col long form:  GeneID <TAB> trans_len <TAB> length

Library size (--libFn; optional, affects RPKM only):
  2 columns:  sampleID <TAB> total_read_count
  Use this to give RPKM the TRUE number of mapped reads. With featureCounts -M,
  the per-column count sum over-counts multi-mapping reads, so the default RPKM
  denominator (column sum) is too large; --libFn overrides it per sample.
  TPM self-normalizes (each column sums to 1e6) and is NOT affected by --libFn.

Output type (--out_type, default tpm):
  tpm | rpkm | both.  With 'both', TPM goes to --outFn/STDOUT and RPKM to --rpkm_out.

Other:
  --skipEleID <file>   extra gene IDs (one per line) to drop. Rows whose ID starts
                       with '__' (htseq/featureCounts summary rows) are always dropped.
  --help

HELP
  exit(1);
}
$o{'help'} and &usage();
defined $o{'cntFn'} or do { warn "[Err] --cntFn is required.\n"; &usage(); };
my $fmt      = lc($o{'fmt'}      // 'auto');
my $out_type = lc($o{'out_type'} // 'tpm');
$out_type =~ m/^(tpm|rpkm|both)$/ or die "[Err] --out_type must be tpm|rpkm|both\n";
if ($out_type eq 'both' and !defined $o{'rpkm_out'}) {
  die "[Err] --out_type both needs --rpkm_out <file> for the RPKM table.\n";
}

# ---- rows to skip ----
my %skip;
if (defined $o{'skipEleID'}) {
  open my $fh, '<', $o{'skipEleID'} or die "[Err] cannot open --skipEleID $o{'skipEleID'}: $!\n";
  while (<$fh>) { chomp; my ($id) = split /\t/, $_; defined $id and $id ne '' and $skip{$id} = 1; }
  close $fh;
}

# ---- read count table ----
open my $cfh, '<', $o{'cntFn'} or die "[Err] cannot open --cntFn $o{'cntFn'}: $!\n";
my $hdr = <$cfh>;
defined $hdr or die "[Err] empty --cntFn\n";
chomp $hdr; $hdr =~ s/\r$//;
my @H = split /\t/, $hdr;

if ($fmt eq 'auto') {
  $fmt = ( @H >= 6 && $H[1] eq 'Chr' && $H[5] eq 'Length' ) ? 'featurecounts' : 'matrix';
}
my ($len_col, $first_cnt_col);   # 0-based
if ($fmt eq 'featurecounts') {
  # locate Length column and take everything after it as samples
  for my $i (0 .. $#H) { $H[$i] eq 'Length' and $len_col = $i; }
  defined $len_col or die "[Err] featurecounts format but no 'Length' column found.\n";
  $first_cnt_col = $len_col + 1;
} elsif ($fmt eq 'matrix') {
  defined $o{'lenFn'} or die "[Err] --fmt matrix requires --lenFn.\n";
  $first_cnt_col = 1;
} else {
  die "[Err] --fmt must be auto|featurecounts|matrix\n";
}
my @samples = @H[$first_cnt_col .. $#H];
@samples or die "[Err] no sample columns found in --cntFn (fmt=$fmt).\n";

# ---- gene lengths ----
my %glen;
if ($fmt eq 'matrix') {
  open my $lf, '<', $o{'lenFn'} or die "[Err] cannot open --lenFn $o{'lenFn'}: $!\n";
  while (<$lf>) {
    m/^\s*(#|$)/ and next;
    chomp; s/\r$//;
    my @a = split /\t/, $_;
    if (@a >= 3 and $a[1] !~ /^[\d.eE+-]+$/) {
      # long form: id, value_type, value ; only trans_len rows are lengths
      $a[1] eq 'trans_len' and $glen{$a[0]} = $a[2];
    } else {
      $glen{$a[0]} = $a[1];
    }
  }
  close $lf;
}

# ---- read counts ----
my (@ids, @len, @cnt);   # @cnt[row] = [c_sample1, c_sample2, ...]
while (<$cfh>) {
  chomp; s/\r$//;
  my @a = split /\t/, $_;
  my $id = $a[0];
  $id =~ /^__/ and next;
  $skip{$id} and next;
  my $L;
  if ($fmt eq 'featurecounts') {
    $L = $a[$len_col];
  } else {
    defined $glen{$id} or die "[Err] no length for gene [$id] in --lenFn.\n";
    $L = $glen{$id};
  }
  ($L and $L > 0) or die "[Err] non-positive length for gene [$id].\n";
  push @ids, $id;
  push @len, $L;
  push @cnt, [ @a[$first_cnt_col .. $#H] ];
}
close $cfh;
@ids or die "[Err] no gene rows read from --cntFn.\n";

# ---- optional true library sizes ----
my %libsize;
if (defined $o{'libFn'}) {
  open my $bf, '<', $o{'libFn'} or die "[Err] cannot open --libFn $o{'libFn'}: $!\n";
  while (<$bf>) {
    m/^\s*(#|$)/ and next;
    chomp; s/\r$//;
    my @a = split /\t/, $_;
    defined $a[1] and $libsize{$a[0]} = $a[1];
  }
  close $bf;
}

# ---- compute RPK, then per-sample sums ----
my $ns = scalar(@samples);
my @rpk;                    # @rpk[row][s]
my @sumRPK = (0) x $ns;     # for TPM denominator
my @colSum = (0) x $ns;     # raw count column sums (default RPKM library size)
for (my $r = 0; $r < @ids; $r++) {
  my $Lkb = $len[$r] / 1000;
  for (my $s = 0; $s < $ns; $s++) {
    my $c = $cnt[$r][$s];
    defined $c and $c ne '' or $c = 0;
    my $v = $c / $Lkb;
    $rpk[$r][$s] = $v;
    $sumRPK[$s] += $v;
    $colSum[$s] += $c;
  }
}

# RPKM library size per sample: --libFn override, else column sum
my @lib;
for (my $s = 0; $s < $ns; $s++) {
  my $sm = $samples[$s];
  $lib[$s] = (defined $libsize{$sm}) ? $libsize{$sm} : $colSum[$s];
  $lib[$s] > 0 or $lib[$s] = 1;   # guard empty samples
  $sumRPK[$s] > 0 or $sumRPK[$s] = 1;
}

# ---- writers ----
sub write_tbl {
  my ($fn, $calc) = @_;   # $calc->($row,$s) -> value
  my $fh;
  if (defined $fn) { open $fh, '>', $fn or die "[Err] cannot write $fn: $!\n"; }
  else             { $fh = \*STDOUT; }
  print {$fh} join("\t", 'GeneID', @samples), "\n";
  for (my $r = 0; $r < @ids; $r++) {
    print {$fh} join("\t", $ids[$r], map { $calc->($r, $_) } 0 .. $ns-1), "\n";
  }
  defined $fn and close $fh;
}
my $tpm_calc  = sub { my ($r,$s)=@_; $rpk[$r][$s] / $sumRPK[$s] * 1e6; };
my $rpkm_calc = sub { my ($r,$s)=@_; $rpk[$r][$s] / ($lib[$s] / 1e6);  };

if    ($out_type eq 'tpm')  { &write_tbl($o{'outFn'}, $tpm_calc); }
elsif ($out_type eq 'rpkm') { &write_tbl($o{'outFn'}, $rpkm_calc); }
else { # both
  &write_tbl($o{'outFn'},   $tpm_calc);
  &write_tbl($o{'rpkm_out'}, $rpkm_calc);
}
