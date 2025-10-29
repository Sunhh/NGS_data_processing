#!/usr/bin/perl
# 20251029: Allow -start_list for gene number
use strict;
use warnings;

-t and !@ARGV and die "perl $0 prefix_of_gene [-start_list in_list.chrN_nextGeneN]  input.gff.JnLoc > map.oriMID_newMID_oriGID_newGID\n";

my $pref=shift;
my %chr2cnt;

if ($ARGV[0] =~ m!^\-start_list$!i) {
  shift;
  my $f1 = shift;
  open F,"<$f1" or die;
  while (<F>) {
    chomp;
    my @ta=split(/\t/, $_);
    $chr2cnt{$ta[0]} = $ta[1];
  }
  close F;
}

my (%onChr);
while (<>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'mrnaID' and next;
  if ($ta[2] =~ m!chr(\d*[1-9]\d*)$!i) {
    push(@{$onChr{$ta[2]}}, [@ta[0,1,6,7,8]]); 
  } elsif ($ta[2] =~ m!chr0+$!i) {
    push(@{$onChr{$ta[2]}}, [@ta[0,1,6,7,8]]); 
  } else {
    die "$_\n";
  }
}

for my $ori_chrID  (sort keys %onChr) {
  $ori_chrID =~ m!chr(\d+)$!i or die;
  my $chrNum = $1; $chrNum += 0;
  @{$onChr{$ori_chrID}} = sort { $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @{$onChr{$ori_chrID}};
  my $cnt=0;
  if (defined $chr2cnt{$chrNum}) {
    $cnt = $chr2cnt{$chrNum}-1;
  }
  for my $t1 (@{$onChr{$ori_chrID}}) {
    $cnt ++;
    my $new_geneID = sprintf("%sC%02dG%06d", $pref, $chrNum, $cnt*10);
    my $new_mrnaID = "$new_geneID.1";
    print join("\t", $t1->[0], $new_mrnaID, $t1->[1], $new_geneID)."\n";
  }
}

# 0	mrnaID	C31_G028498
# 1	geneID	C31_G028498-G
# 2	SeqID	C31_Chr16
# 3	mrnaStart	23350
# 4	mrnaEnd	24470
# 5	mrnaStrand	-
# 6	CDSStart	23350
# 7	CDSEnd	24470
# 8	LenInCDS	681
# 9	CDSBlocks	24096,24470;23951,23997;23509,23697;23350,23419
# 10	CDSBlocksNum	4
# 11	5UTR	0
# 12	3UTR	0
