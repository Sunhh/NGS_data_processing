#!/usr/bin/perl
# 261001: It seems that GenomeScope2.0 has changed the estimates in PDF, so they are not "max" in summary.txt any longer.
# v2 (261001, Claude Code helped from extract_genomescope_summary.pl): report the values shown in the GenomeScope plots.
#   summary.txt gives min/max = estimate -/+ 2 SE; which values the plots show depends on the GenomeScope version.
#   genome_size, unique_region, repeat_region:
#     GenomeScope 1.0, and 2.0 run before 2023-11-23 (commit 0faf3f8, e.g. release v1.0.0): the max column (kcov - 2 SE).
#     GenomeScope 2.0 run from 2023-11-23 (e.g. release v2.0.1): the estimate. A length is k-mers / (p * kcov),
#       so the estimate (kcov = kmercov) is the harmonic mean of min and max.
#     Both 2.0 versions write "GenomeScope version 2.0" in summary.txt; give -gs2_len max for the older one.
#   hete%: the estimate in every version.
#     GenomeScope 1.0: midpoint of Heterozygosity min/max (the old script reported the max).
#     GenomeScope 2.0: p=2-4, midpoint of Heterozygous min/max; if a bound was clipped at 0% or 100%, the sum of the
#       r estimates in model.txt. p>=5, the single value in summary.txt (already the estimate). p=1, N/A.
#   kmer_depth: p * kmercov in model.txt (2 * kmercov for 1.0 and for p=2, as before).
#   model.txt is the one next to summary.txt with the same name prefix (BC19_summary.txt -> BC19_model.txt).
use strict;
use warnings;
use fileSunhh;
use Getopt::Long;
use List::Util qw(sum);

my %opts = ('gs2_len' => 'est');
GetOptions(\%opts, "gs2_len=s", "help!") or die "[Err] Bad options. See: perl $0 -help\n";
my $usage = <<"HH";
perl $0 [-gs2_len est|max] out.of.genomescope/BC19/summary.txt out.of.genomescope/BC20/summary.txt > summary_gscope.tbl

  -gs2_len est   [default] GenomeScope 2.0 run from 2023-11-23 (release v2.0.1 or later): lengths are the estimates.
  -gs2_len max   GenomeScope 2.0 run before 2023-11-23 (e.g. release v1.0.0): lengths are the max column, as in 1.0.
  GenomeScope 1.0 is recognised from the first line of summary.txt and does not need -gs2_len.
HH
($opts{'help'} or !@ARGV) and die $usage;
$opts{'gs2_len'} =~ m!^(est|max)$! or die "[Err] -gs2_len should be est or max\n\n$usage";

my $num = qr/[-+]?(?:\d[\d,]*(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?/; # numbers as GenomeScope writes them: with "," or in scientific notation
my @outKey = qw(filename k k_depth genome_hap hetePerc genome_uni genome_rep);
print STDOUT join("\t", qw/filename kmer kmer_depth genome_size hete% unique_region repeat_region/)."\n";
for my $fn (@ARGV) {
  open F,'<',"$fn" or die "[Err] Failed to open file $fn\n";
  my %h;
  $h{'filename'} = $fn;
  my $version = '';
  my $p = 2;      # GenomeScope 1.0 summary.txt has no "p = " line
  my %len;        # Haploid/Repeat/Unique => [min, max] in bp
  my @het;        # [min, max] of Heterozygosity / Heterozygous (...)
  my $hetSingle;  # p>=5: one value
  my @hetClass;   # p=3,4: [min, max] of aab, abc, ...
  while (<F>) {
    chomp;
    if (m!^GenomeScope version (\S+)!) {
      $version = $1;
    } elsif (m!^p\s*=\s*(\d+)\s*$!) {
      $p = $1;
    } elsif (m!^k\s*=\s*(\d+)\s*$!) {
      $h{'k'} = $1;
    } elsif (m!^Genome (Haploid|Repeat|Unique) Length\s+($num)\s*bp\s+($num)\s*bp!) {
      $len{$1} = [ &_n($2), &_n($3) ];
    } elsif (m!^(?:Heterozygosity|Heterozygous\s*\([^)]+\))\s+($num)\%\s+($num)\%\s*$!) {
      @het = ( &_n($1), &_n($2) );
    } elsif (m!^(?:Heterozygosity|Heterozygous\s*\([^)]+\))\s+($num)\%\s*$!) {
      $hetSingle = &_n($1);
    } elsif (m!^[a-f]+\s+($num)\%\s+($num)\%\s*$!) {
      push @hetClass, [ &_n($1), &_n($2) ];
    }
  }
  close F;
  my $isV1 = ($version =~ m!^1\.!);
  $isV1 or $version =~ m!^2\.! or warn "[Wrn] $fn: ".($version eq '' ? "no 'GenomeScope version' line" : "unknown GenomeScope version $version").", read as 2.0\n";

  my $absFn   = &fileSunhh::_abs_path($fn);
  my $dirName = &fileSunhh::_dirname( $absFn );
  my $modelFn = $absFn;
  $modelFn =~ s!summary\.txt$!model.txt! or $modelFn = "$dirName/model.txt";
  my %est; # estimates in model.txt: kmercov, r (GenomeScope 1.0) or r1, r2, ... (2.0)
  if (-e $modelFn) {
    open F2,'<',"$modelFn" or die "[Err] Failed to open file $modelFn\n";
    while (<F2>) {
      m!^(kmercov|r\d*)\s+($num)\s! and $est{$1} = $2;
    }
    close F2;
  }

  ## min and max are -1 bp when the model failed to converge
  if (defined $len{'Haploid'} and $len{'Haploid'}[0] > 0 and $len{'Haploid'}[1] > 0) {
    my %lenKey = qw(Haploid genome_hap Repeat genome_rep Unique genome_uni);
    for my $t (keys %lenKey) {
      defined $len{$t} or next;
      my $v = ($isV1 or $opts{'gs2_len'} eq 'max') ? $len{$t}[1] : &_harmonic_mean( @{$len{$t}} );
      defined $v and $v >= 0 and $h{$lenKey{$t}} = sprintf("%.0f", $v);
    }
    if ($isV1) {
      @het and $h{'hetePerc'} = sprintf("%.6g", ($het[0]+$het[1])/2);
    } elsif ($p >= 5) {
      defined $hetSingle and $h{'hetePerc'} = sprintf("%.6g", $hetSingle);
    } elsif ($p >= 2 and @het) {
      my @comp = (@hetClass) ? @hetClass : ([@het]);
      if (grep { $_->[0] <= 0 or $_->[1] >= 100 } @comp) {
        my @r = grep { m!^r\d+$! } keys %est;
        if (@r) {
          $h{'hetePerc'} = sprintf("%.6g", 100 * sum(@est{@r}));
          warn "[Wrn] $fn: a Heterozygous bound is clipped at 0% or 100%, so hete% is the sum of r estimates in $modelFn\n";
        } else {
          warn "[Wrn] $fn: a Heterozygous bound is clipped at 0% or 100% and $modelFn has no r estimates, hete% is N/A\n";
        }
      } else {
        $h{'hetePerc'} = sprintf("%.6g", ($het[0]+$het[1])/2);
      }
    }
    defined $est{'kmercov'} and $h{'k_depth'} = $est{'kmercov'} * $p;
  }
  for (@outKey) {
    $h{$_} //= "N/A";
  }
  print STDOUT join("\t", @h{@outKey})."\n";
}

sub _n {
  my $v = shift;
  $v =~ s!,!!g;
  return $v + 0;
}

## Both bounds are k-mers / (p * kcov) with kcov = kmercov +/- 2 SE, so the estimate (kcov = kmercov) is 2/(1/min + 1/max).
sub _harmonic_mean {
  my ($lo, $hi) = @_;
  ($lo == 0 and $hi == 0) and return 0;
  ($lo > 0 and $hi > 0) or return undef;
  return 2 / (1/$lo + 1/$hi);
}

# Sunhh@swift:/data/Sunhh/wmhifi/genome_size$ less -S out.of.genomescope/BC19/summary.txt 
# GenomeScope version 1.0
# k = 81
# 
# property                      min               max               
# Heterozygosity                0.0523306%        0.0526906%        
# Genome Haploid Length         414,606,094 bp    414,694,212 bp    
# Genome Repeat Length          79,935,019 bp     79,952,008 bp     
# Genome Unique Length          334,671,075 bp    334,742,204 bp    
# Model Fit                     97.7185%          98.7888%          
# Read Error Rate               0.226308%         0.226308%         
