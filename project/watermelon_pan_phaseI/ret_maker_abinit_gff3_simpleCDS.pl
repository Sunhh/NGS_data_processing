#!/usr/bin/perl
# [2/25/2022] I'd like to add _QI parse to correctly divide 'match_part' proportion in the future. Currently I just deal with simple cases.
use strict;
use warnings;
use fileSunhh;

!@ARGV and die "perl $0 ab_init.ID_list maker_all.gff3 > ab_init.match.gff3\n";

my $f1_list = shift;
my $f2_gff3 = shift;

my %h;
for (&fileSunhh::load_tabFile($f1_list)) {
  $_->[0] =~ s!\-processed\-gene\-!-abinit-gene-!;
  $h{$_->[0]} = 1;
}
my $ofh2 = &openFH($f2_gff3, '<');
while (<$ofh2>) {
  chomp;
  m!^\s*#! and next;
  m!^\>! and last;
  my @ta=split(/\t/, $_);
  my $is_out = 0;
  # $ta[1] =~ m!^(repeatmasker|protein_gff:blastx)$!i and next;
  if ($ta[2] =~ m!^match$!i) {
    $ta[8] =~ s!^ID=[^\s]+;Name=([^\s;]+)!ID=$1!i or die "[Err] Bad ta8 [$ta[8]] in line: [$_]\n";
    if (defined $h{$1}) {
      my $mID = $1;
      my $gID = $mID; $gID =~ s!\-mRNA\-\d+$!! or die "[Err] gID [$gID]\n";
      $ta[8] =~ m!_QI=0\|[^\s]+\|0\|\d+;! or die "[Err] Bad _QI in ta8 [$ta[8]] in line [$_]\n";
      print join("\t", @ta[0,1], "gene", @ta[3,4], ".", @ta[6,7], "ID=$gID")."\n";
      print join("\t", @ta[0,1], "mRNA", @ta[3..8])."\n";
    }
  } elsif ($ta[2] =~ m!^match_part$!i) {
    $ta[8] =~ s!^(ID=[^\s]+;)Target=([^\s;]+)!${1}Parent=$2;Target=$2!i or die "[Err] Bad ta8 [$ta[8]] in line: [$_]\n";
    if (defined $h{$2}) {
      my $mID = $2;
      print join("\t", @ta[0,1], "exon", @ta[3..7], "Parent=$mID")."\n";
      print join("\t", @ta[0,1], "CDS", @ta[3..8])."\n";
    }
  }
}
close($ofh2);

# [Sunhh@panda iprscan]$ more test.gff3 
# s030831 augustus_masked match   88      900     0.41    +       .       ID=s030831:hit:20299:4.5.0.0;Name=augustus_masked-s030831-abinit-gene-0.3-mRNA-1;_AED=1.00;_eAED=1.00;_QI=0|0|0|0|1|1|3|0|74;target_length=1356
# s030831 augustus_masked match_part      88      152     0.41    +       .       ID=s030831:hsp:43060:4.5.0.0;Parent=s030831:hit:20299:4.5.0.0;Target=augustus_masked-s030831-abinit-gene-0.3-mRNA-1 1 65 +;Gap=M65
# s030831 augustus_masked match_part      511     603     0.41    +       .       ID=s030831:hsp:43061:4.5.0.0;Parent=s030831:hit:20299:4.5.0.0;Target=augustus_masked-s030831-abinit-gene-0.3-mRNA-1 66 158 +;Gap=M93
# s030831 augustus_masked match_part      834     900     0.41    +       .       ID=s030831:hsp:43062:4.5.0.0;Parent=s030831:hit:20299:4.5.0.0;Target=augustus_masked-s030831-abinit-gene-0.3-mRNA-1 159 225 +;Gap=M67

