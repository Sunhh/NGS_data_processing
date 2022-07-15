#!/usr/bin/perl
## https://en.wikipedia.org/wiki/BED_(file_format)#:~:text=The%20BED%20(Browser%20Extensible%20Data,adopted%20by%20other%20sequencing%20projects.
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "help!",
  "for_nucdiff!",
);

-t and !@ARGV and die "perl $0 -for_nucdiff ref_struct.gff > ref_struct.bed\n";

while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  my $id="";
  if ($opts{'for_nucdiff'}) {
    $ta[8] =~ m!(?:^|\s|;)ID=([^\s;]+)!i and $id=$1;
    my $type="";
    $ta[8] =~ m!(?:^|\s|;)Name=([^\s;]+)!i and $type=$1;
    my $xBlkLen=-1;
    if ($ta[8] =~ m!(?:^|\s|;)(?:overlap_len|subst_len|del_len|ins_len|repeat_len|region_len|duplic_len|blk_length|blk_len)=([\-+\d.]+)!) {
      $xBlkLen = $1;
    } elsif ( $type =~ m!relocation|translocation! ) {
      # $xBlkLen = 0;
    } else {
      die "Failed to get length in formation:$_\n";
    }
    # $id="ID-${id}:Type-$type:xBlkLen-$xBlkLen";
    $id="${id}:$type:$xBlkLen";
  }
  if ($id eq "") {
    print join("\t", $ta[0], $ta[3]-1, $ta[4])."\n";
  } else {
    print join("\t", $ta[0], $ta[3]-1, $ta[4], $id)."\n";
  }
}

# ##gff-version 3
# ##sequence-region	C39_Chr01	1	26126445
# C39_Chr01	NucDiff_v2.0	SO:0001873	1	1	.	.	.	ID=SV_1.2;Name=translocation-overlap;overlap_len=4886;query_sequence=C31_Chr01;query_coord=9464;breakpoint_query=9464-14349;blk_query=9464-11510801;blk_ref=1-12001526;blk_query_len=11501338;blk_ref_len=12001526;color=#A0A0A0
# C39_Chr01	NucDiff_v2.0	SO:1000002	174594	175294	.	.	.	ID=SV_2;Name=substitution;subst_len=701;query_dir=1;query_sequence=C31_Chr01;query_coord=184573-185273;color=#42C042
# Number	KeyCols
# 678	overlap_len
# 5327	subst_len
# 3478	del_len
# 8107	ins_len
# 385	blk_len
# 22	blk_query_len
# overlap_len
# repeat_len
# region_len
# duplic_len
