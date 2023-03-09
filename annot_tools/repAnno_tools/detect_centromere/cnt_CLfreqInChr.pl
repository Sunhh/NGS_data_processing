#!/usr/bin/perl
use strict;
use warnings;

my $htxt = <<HHH;

perl $0 \\
 px_pyTanFinder \\
 list_chrID.long_short \\
 out_pyTanFinder/out_merger/merger_output_XXX/after_clusterin_summary_XXX.txt \\
> out_pyTanFinder/out_merger/merger_output_XXX/after_clusterin_summary_XXX.txt.cntChr

list_chrID.long_short : chrID_in_seqfile  \\t  chrID_like_Chr01 \\n


HHH

-t and scalar(@ARGV) < 3 and die "$htxt";
my $min_monomer_len = 70; # Set as 70 because I think most centromeric repeat unit should size 80-200 bp.
my $px = shift;
my $f1List = shift;
open F1,'<',"$f1List" or die;
my (%idL2S);
while (<F1>) {
  m!^\s*(#|$)! and next;
  chomp;
  my @ta=split(/\t/, $_);
  $idL2S{$ta[0]} = $ta[1];
}
close F1;

my @trs;
while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  $ta[0] eq 'CL-' and next;
  $ta[0] eq 'Cluster' and next;
  $ta[3] = &cnvt2number($ta[3]);
  $ta[4] = &cnvt2number($ta[4]);
  $ta[3] < $min_monomer_len and next;
  # $ta[0] eq "CL360" and die "$_\n$ta[3]\n";
  my %chrs;
  for my $a1 (split(/,/, $ta[6])) {
    my @tb = split(/\<\*\>/, $a1);
    my ($chrID, $repN, $monoLen) = @tb[0,2,3];
    $chrID =~ s!$px!!; defined $idL2S{$chrID} or next;
    my $chrIDS = $idL2S{$chrID};
    # $a1 =~ m!^$px${chrID}\<\*\>\d+\<\*\>([\d.]+)\<\*\>(\d+)$! or die "bad ID [$a1]\n";
    $monoLen < $min_monomer_len and next;
    if (defined $chrs{$chrIDS}) {
      $chrs{$chrIDS}[0] < $repN and $chrs{$chrIDS} = [ $repN, $a1 ];
    } else {
      $chrs{$chrIDS} = [ $repN, $a1 ];
    }
  }
  my @id_chr = sort keys %chrs;
  my $n_chr = scalar(@id_chr);
  push(@trs, [$ta[0], $ta[1], $ta[3], $ta[4], $n_chr, join(";", @id_chr), join(";", map { $_->[1] } @chrs{@id_chr})]);
}
print STDOUT join("\t", qw/CL_ID Repres_ID monomer_len accumu_repeat_times cover_chr_num cover_chr_IDs repres_in_chrs/)."\n";
for my $ar1 (sort { $b->[3]<=>$a->[3] || $b->[4]<=>$a->[4] } @trs) {
  print STDOUT join("\t", @$ar1)."\n";
}

sub cnvt2number {
  my ($txt) = @_;
  if ($txt =~ m!^\d+$!) {
    return($txt + 0);
  } elsif ($txt =~ m!^(\d+)\,(\d+)$!) {
    my $v1 = $1;
    my $v2 = "0.".$2;
    return($v1 + $v2);
  } elsif ($txt =~ m!^(\d+)\.(\d+)$!) {
    return($txt + 0);
  } else {
    die "[Err] Wrong number [$txt]\n";
  }
}

# 0       Cluster                 CL1     CL2
# 1       ID                      XXXScaffold7014<*>8518236<*>6.4<*>201   XXXChr6.1<*>4434616<*>5.5<*>146
# 2       Sequence                CACGAGCTTCCGCTCATGCCTTCCCGGCAACCCTTGAGCTTCCGCTCACGAGCTTCCCGCTCATGCCTTCCCGGCACCCCGAGCTTC>
# 3       Monomer length          201,0   146,0
# 4       Accumulated abundance   1876922,5000000102      98991,10000000003
# 5       Number of connections   3196    253
# 6         xx                    XXXChr6.1<*>4496672<*>5.1<*>25,XXXChr7.4<*>6656324<*>8.0<*>26,XXXChr3.4<*>1713650<*>7.2<*>26,XXXScaffold39<*>7116737<*>8.0<*>26,XXXChr7.2<*>5734251<*>23.2<*>26,XXXChr4.2<*>3470337<*>5.8<>


