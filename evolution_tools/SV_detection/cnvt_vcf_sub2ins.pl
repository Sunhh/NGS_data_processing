#!/usr/bin/perl
# 6/6/2023: Change long SUBstitutions to INSertions on the left. I don't include DELetion because it cannot provide true junction sites and can introduce duplicated reference positions/alleles.
use strict;
use warnings;

!@ARGV and die "perl $0 maxLen_5 in.vcf > fixed.vcf\n";

my $maxLen = shift;

while (<>) {
  chomp;
  m!^#! and do { print "$_\n"; next; };
  my @ta=split(/\t/, $_);
  $ta[4] eq '<INV>' and do { print "$_\n"; next; };
  if ( (length($ta[3]) > $maxLen and length($ta[4]) > 1) or (length($ta[3]) > 1 and length($ta[4]) > $maxLen) ) {
    $ta[3] = substr($ta[3], 0, 1);
    $ta[4] = $ta[3] . $ta[4];
    print join("\t", @ta)."\n";
  } else {
    print "$_\n";
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


