#!/usr/bin/perl
# 2/1/2023: Count deletion/insertion in vcf file. Always use the first allele to determine SV type.
#           Allele label '.' is not accepted.
# 3/9/2023: Update the code, providing same results. Sizes for -small_indel_max filtering:
#             alternative allele: longest - alternative allele;
#             MNP : length(ref_allele)
#             SV size 1 : index(base_ref, base_alt)==0: len_ref-len_alt;
#             SV size 2 : index(base_alt, base_ref)==0: len_alt-len_ref;
#             SV size 3 : max(len_ref, len_alt)
# 6/7/2023: Skip <INV> allele.

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "class_limit:s", # '1-19 20-49 50-999  1000-9999 10000-49999 50000-99999 100000-199999 200000-299999 300000-'
  "small_indel_max:i", # 0.
"help!");

sub usage {
  print <<HH;
################################################################################
bcftools view 97103v2.aw_ref_snps.vcf.gz | perl $0 > 97103v2.aw_ref_snps.vcf.gz.sv_summary\n

  -class_limit         '1-19 20-49 50-999  1000-9999 10000-49999 50000-99999 100000-199999 200000-299999 300000-'
  -small_indel_max     [0]

############
Always use the first allele to determine SV type.
  Insertion (INS): length(ref) > length(alt);
    case 1: A=>AAT: INS_size=2;
    case 2: T=>AAT: INS_size=3; DEL_size=1;
  Deletion  (DEL): length(ref) < length(alt);
    case 1: AAT=>AA: DEL_size=1;
    case 2: AAT=>TT: DEL_size=3; INS_size=2;
  Substitution (SUB): length(ref) == length(alt); SUB=INS+DEL;
################################################################################
HH
  exit(1);
}

-t and !@ARGV and &usage();
$opts{'help'} and &usage();
$opts{'small_indel_max'} //= 0;

# Report information:
### Insertion(INS)/Deletion(DEL)/Substitution(SUB):
###   Total number and length;
###   Summary in each class limits.
###   SUB = INS + DEL

my $cls_lim_txt = '1-19 20-49 50-999  1000-9999 10000-49999 50000-99999 100000-199999 200000-299999 300000-';
defined $opts{'class_limit'} and $opts{'class_limit'} ne '' and $cls_lim_txt = $opts{'class_limit'};
my $max_cls_lim_txt_len = -1;
my @cls_lim_arr;
for my $a1 (split(/\s+/, $cls_lim_txt)) {
  $max_cls_lim_txt_len < length($a1) and $max_cls_lim_txt_len = length($a1);
  if      ($a1 =~ m!^(\d+)\-(\d+)$!) {
    push(@cls_lim_arr, [$1, $2, $a1]);
  } elsif ($a1 =~ m!^(\d+)\-$!) {
    push(@cls_lim_arr, [$1, -1, $a1]);
  } else {
    die "failed to parse |$a1|\n";
  }
}
@cls_lim_arr = sort { $a->[0] <=> $b->[0] } @cls_lim_arr;
$max_cls_lim_txt_len ++;


my %cls_cnt;
while (<>) {
  m!^\s*#|^\s*$! and next;
  chomp;
  my @ta=split(/\t/, $_);
  $ta[3] = uc($ta[3]);
  $ta[4] = uc($ta[4]);
  $ta[3] =~ m!^[ATGCN]+$! or die "[Err] allele should be ATGCN. here is |$ta[3]|\n";
  my $base_ref = $ta[3];
  my $len_ref = length($ta[3]);
  # Fix to use the first allele.
  my @tb = split(/,/, $ta[4]);
  my ($base_alt, $len_alt);
  for my $b1 (@tb) {
    $b1 eq '<INV>' and next;
    $base_alt //= $b1;
    $len_alt //= length($b1);
    $len_alt < length($b1) and do { $len_alt=length($b1); $base_alt=$b1; };
  }
  defined $base_alt or next;
  $base_ref = uc($base_ref);
  $base_alt = uc($base_alt);
  unless ( $base_alt =~ m!^[ATGCN]+$!i) {
    warn "[Err] allele should be ATGCN. here is |$base_alt! in: @tb\n";
    next;
  }

  if ($len_ref == 1 and $len_alt == 1) {
    &add_snp(1);
    next;
  }
  if      (index($base_ref, $base_alt) == 0) {
    if ($len_ref-$len_alt <= $opts{'small_indel_max'}) {
      &add_small_del($len_ref-$len_alt);
    } else {
      &add_del($len_ref-$len_alt);
    }
    if ($ta[7] =~ m!(SVTYPE=[^;]+)! and $ta[7] !~ m!SVTYPE=(?:DEL|SUB|CPX)\s*(?:;|$)!) {
      warn "[Wrn] Expected DEL/SUB SV at [$ta[0] $ta[1]] [refLen=$len_ref altLen=$len_alt] $1\n";
    }
  } elsif (index($base_alt, $base_ref) == 0) {
    if ($len_alt-$len_ref <= $opts{'small_indel_max'}) {
      &add_small_ins($len_alt-$len_ref);
    } else {
      &add_ins($len_alt-$len_ref);
    }
    if ($ta[7] =~ m!(SVTYPE=[^;]+)! and $ta[7] !~ m!SVTYPE=(?:INS|SUB|CPX)\s*(?:;|$)!) {
      warn "[Wrn] Expected INS/SUB SV at [$ta[0] $ta[1]] [refLen=$len_ref altLen=$len_alt] $1\n";
    }
  } elsif ($len_ref <= $opts{'small_indel_max'} and $len_alt <= $opts{'small_indel_max'}) {
    &add_small_del($len_ref);
    &add_small_ins($len_alt);
  }else {
    &add_del($len_ref);
    &add_ins($len_alt);
  }
}

print STDOUT join("\t", qw/Class_limit  INS_cnt INS_len DEL_cnt DEL_len SNP_cnt SNP_len INSsmall_cnt INSsmall_len DELsmall_cnt DELsmall_len/)."\n";
for my $c1 ((map {$_->[2]} @cls_lim_arr), 'all') {
  my @o1 = (sprintf("%-${max_cls_lim_txt_len}s", $c1));
  for my $type1 (qw/ins del snp small_ins small_del/) {
    $cls_cnt{$type1}{$c1} //= { 'num'=> 0, 'len' =>0 };
    push(@o1, $cls_cnt{$type1}{$c1}{'num'}, $cls_cnt{$type1}{$c1}{'len'});
  }
  print STDOUT join("\t", @o1)."\n";
}


sub add_snp {
  my ($size) = @_;
  for my $c1 (@cls_lim_arr) {
    $c1->[0] <= $size or last;
    $c1->[1] >= $size or $c1->[1] == -1 or next;
    $cls_cnt{'snp'}{$c1->[2]}{'num'} ++;
    $cls_cnt{'snp'}{$c1->[2]}{'len'} += $size;
  }
  $cls_cnt{'snp'}{'all'}{'num'} ++;
  $cls_cnt{'snp'}{'all'}{'len'} ++;
  return();
}# add_ins()

sub add_ins {
  my ($size) = @_;
  for my $c1 (@cls_lim_arr) {
    $c1->[0] <= $size or last;
    $c1->[1] >= $size or $c1->[1] == -1 or next;
    $cls_cnt{'ins'}{$c1->[2]}{'num'} ++;
    $cls_cnt{'ins'}{$c1->[2]}{'len'} += $size;
  }
  $cls_cnt{'ins'}{'all'}{'num'} ++;
  $cls_cnt{'ins'}{'all'}{'len'} += $size;
  return();
}# add_ins()

sub add_small_ins {
  my ($size) = @_;
  for my $c1 (@cls_lim_arr) {
    $c1->[0] <= $size or last;
    $c1->[1] >= $size or $c1->[1] == -1 or next;
    $cls_cnt{'small_ins'}{$c1->[2]}{'num'} ++;
    $cls_cnt{'small_ins'}{$c1->[2]}{'len'} += $size;
  }
  $cls_cnt{'small_ins'}{'all'}{'num'} ++;
  $cls_cnt{'small_ins'}{'all'}{'len'} += $size;
  return();
}# add_small_ins()

sub add_small_del {
  my ($size) = @_;
  for my $c1 (@cls_lim_arr) {
    $c1->[0] <= $size or last;
    $c1->[1] >= $size or $c1->[1] == -1 or next;
    $cls_cnt{'small_del'}{$c1->[2]}{'num'} ++;
    $cls_cnt{'small_del'}{$c1->[2]}{'len'} += $size;
  }
  $cls_cnt{'small_del'}{'all'}{'num'} ++;
  $cls_cnt{'small_del'}{'all'}{'len'} += $size;
  return();
}# add_small_del()

sub add_del {
  my ($size) = @_;
  for my $c1 (@cls_lim_arr) {
    $c1->[0] <= $size or last;
    $c1->[1] >= $size or $c1->[1] == -1 or next;
    $cls_cnt{'del'}{$c1->[2]}{'num'} ++;
    $cls_cnt{'del'}{$c1->[2]}{'len'} += $size;
  }
  $cls_cnt{'del'}{'all'}{'num'} ++;
  $cls_cnt{'del'}{'all'}{'len'} += $size;
  return();
}# add_del()


