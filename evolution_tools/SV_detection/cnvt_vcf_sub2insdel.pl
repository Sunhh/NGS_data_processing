#!/usr/bin/perl
# 6/6/2023: Change long SUBstitutions to INSertions on the left. I don't include DELetion because it cannot provide true junction sites and can introduce duplicated reference positions/alleles.
# 8/9/2023: Update SVLEN and END. In the past, we determined to record only large insertions from large substitutions, but now I want to record large deletions, because without these deletions, the sizes of inserted bp and deleted bp are imbalanced and this problem exists after switching Ref and Qry genomes. 'SVTYPE=SUB' kept unchanged.
# 8/11/2023: Fix a bug which duplicates a variant of original INS/DEL.
# 8/17/2023: Correct '.' alleles in REF and ALT.
use strict;
use warnings;
use fastaSunhh;
my $fs_obj = fastaSunhh->new();

!@ARGV and die "perl $0 r.fa maxLen_5 in.vcf > fixed.vcf\n";

my $refFaFn = shift;
my $maxLen = shift;

my %refSeq = %{ $fs_obj->save_seq_to_hash( 'faFile'=>$refFaFn ) };
for my $k (keys %refSeq) { $refSeq{$k}{'seq'} =~ s!\s!!g; $refSeq{$k}{'seq'} = uc($refSeq{$k}{'seq'}); }

while (<>) {
  chomp;
  m!^#! and do { print STDOUT "$_\n"; next; };
  my @ta=split(/\t/, $_);
  $ta[4] eq '<INV>' and do { print STDOUT "$_\n"; next; };
  if ($ta[3] eq '.') {
    $ta[3] = substr($refSeq{$ta[0]}{'seq'}, $ta[1]-1, 1);
    $ta[4] = $ta[3] . $ta[4];
    print STDOUT join("\t", @ta)."\n";
    next;
  }
  if ($ta[4] eq '.') {
    $ta[1] --;
    $ta[4] = substr($refSeq{$ta[0]}{'seq'}, $ta[1]-1, 1);
    $ta[3] = $ta[4] . $ta[3];
    print STDOUT join("\t", @ta)."\n";
    next;
  }
  my $del_len = length($ta[3]);
  my $ins_len = length($ta[4]);
  $del_len == 1 and $ins_len == 1 and do { print "$_\n"; next; };
  if (($del_len > $maxLen and $ins_len > 1) or ($ins_len > $maxLen and $del_len > 1)) {
    # Deleted sequence.
    my @tb = @ta;
    $tb[1] --;
    $tb[4] = substr($refSeq{$tb[0]}{'seq'}, $tb[1]-1, 1);
    $tb[3] = $tb[4] . $tb[3];
    my @tc = split(/;/, $tb[7]);
    for my $a1 (@tc) {
      $a1 =~ s!^SVLEN=([+-]?\d+)$!SVLEN=-$del_len!;
    }
    $tb[7] = join(";", @tc);
    $tb[7] =~ m!SVTYPE=[^\s;]! or $tb[7] .= ";SVTYPE=DEL:MNP";
    print STDOUT join("\t", @tb)."\n";

    # Inserted sequence.
    my @td = @ta;
    ### Update REF and ALT.
    $td[1] --;
    $td[3] = substr($refSeq{$td[0]}{'seq'}, $td[1]-1, 1);
    $td[4] = $td[3] . $td[4];
    ### Update SVLEN and END;
    my @te = split(/;/, $td[7]);
    # ALGORITHMS=NucDiff_struct;SVTYPE=SUB;END=23085209;SVLEN=-4747
    my @tf;
    for my $a1 (@te) {
      $a1 =~ m!^END=\d+! and next;
      $a1 =~ s!^SVLEN=([+-]?\d+)$!SVLEN=$ins_len!;
      push(@tf, $a1);
    }
    $td[7] = join(";", @tf);
    $td[7] =~ m!SVTYPE=[^\s;]! or $td[7] .= ";SVTYPE=INS:MNP";
    print STDOUT join("\t", @td)."\n";
  } else {
    print STDOUT "$_\n";
  }
}

# 22CEXU43_Chr01  171772  .       CA      AG      .       .       ALGORITHMS=NucDiff_snps GT      1/1
# 22CEXU43_Chr01  22540927        .       T       <INV>   .       PASS    ALGORITHMS=NucDiff_struct;SVTYPE=INV;END=22542529;SVLEN=
#
# Take 22CEXU6 as an example:
#
### It should be OK to accept 
### Of the 7136123 ref_snps, 
###   168 variants have their length(REF) or length(ALT) > 10 bp and the other length > 1 bp.
###   334371 variants have both REF and ALT larger than 1 bp.
###     1784 variants have one > 5 bp and the other > 1 bp.
##### awk ' $5 != "<INV>" && ((length($4) > 10 && length($5) >1 ) || (length($4) > 1 && length($5) > 10 )) '
##### awk ' $5 != "<INV>" && ((length($4) > 1 && length($5) >1 ) || (length($4) > 1 && length($5) > 1 )) '
##### awk ' $5 != "<INV>" && ((length($4) > 5 && length($5) >1 ) || (length($4) > 1 && length($5) > 5 )) '
### Of the 3959 ref_struct, 
###   1292 variants have one > 10 bp and the other > 1 bp.
###   1307 variants have one >  5 bp and the other > 1 bp.
### I prefer to change all these long (> 5 bp) SUBstitutions to INSertions after the left REF boundary, withouting changing the 'SVTYPE'. In this way, it is easier to include sequences that are absent from the reference and may help shorten the mapping time.


