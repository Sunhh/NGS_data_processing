#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use fileSunhh;
use fastaSunhh;
my $fs=fastaSunhh->new();
use LogInforSunhh;
my %opts;
GetOptions(\%opts,
  "help!",
  "fn_grpList:s",
);

my $help_txt = <<HH;
perl $0 OG_aln.sep_cds.fa > OG_aln.1_cds.fa

 -fn_grpList    [] 'C01_1cp.HOG_grp_list';

# OG_aln.sep_cds.fa is the output of script 02.list_run_muscle.pl;

HH

!@ARGV and &LogInforSunhh::usage($help_txt);

my $inFas = shift;
my %s1=%{$fs->save_seq_to_hash("faFile"=>$inFas)};
my (@tax, %alnID_inTax, %geneID2alnID);
for my $k1 (sort { $s1{$a}{'Order'} <=> $s1{$b}{'Order'} } keys %s1) {
  $s1{$k1}{'head'} =~ m!^(\S+)\.(\d+)\s(?:\S+\s+)*\[(\S+)\](\s\d+)?\s*$! or die "[Err] k1 [$k1] header: $s1{$k1}{'head'}\n";
  my ($taxID, $gRank, $geneID) = ($1, $2, $3);
  $s1{$k1}{'seq'} =~ s!\s!!g;
  defined $alnID_inTax{$taxID} or push(@tax, $taxID);
  push(@{$alnID_inTax{$taxID}}, [$taxID, $gRank]);
  defined $geneID2alnID{$taxID}{$geneID} and die "[Err] I don't accept repeated gene IDs now. [$geneID]\n";
  $geneID2alnID{$taxID}{$geneID} = "$taxID.$gRank";
}

my $has_grp = 0;
my @tax2gene = ();
my @tax2txt  = ();

if (defined $opts{'fn_grpList'}) {
  my @taxID_txt;
  for my $a1 ( &fileSunhh::load_tabFile($opts{'fn_grpList'}) ) {
    for (my $i=0; $i<@$a1; $i++) {
      push(@{$tax2gene[$i]}, $a1->[$i]);
      for my $t1 (@tax) {
        defined $geneID2alnID{$t1}{$a1->[$i]} or next;
        $taxID_txt[$i]{$t1} ++;
        last;
      }
    }
  }
  for my $a2 ( @taxID_txt ) {
    if (!(defined $a2) or scalar(keys %$a2) == 0) {
      $a2 = ""; next;
    }
    ($a2) = sort { $a2->{$b} <=> $a2->{$a} } keys %$a2;
  }
  my (@a3, @a4);
  for (my $i=0; $i<@tax2gene; $i++) {
    $taxID_txt[$i] eq '' and next;
    #for my $g1 (@{$tax2gene[$i]}) {
    #  defined $geneID2alnID{$taxID_txt[$i]}{$g1} or die "[Err] No alnID found for [$g1]\n";
    #}
    # push(@a3, [map { $geneID2alnID{$taxID_txt[$i]}{$_} } @{$tax2gene[$i]}]);
    push(@a3, [ @{$tax2gene[$i]} ]);
    push(@a4, $taxID_txt[$i]);
  }
  # @tax2gene  = @a3;
  @tax2txt = @a4;
  @tax2gene = ();
  for (my $j=0; $j<@{$a3[0]}; $j++) {
    my $is_allAln = 1;
    my @toadd;
    for (my $i=0; $i<@tax2txt; $i++) {
      # defined $geneID2alnID{$tax2txt[$i]}{$a3[$i][$j]} or die "[Err] No alnID found for [{$tax2txt[$i]}{$a3[$i][$j]}]\n";
      defined $geneID2alnID{$tax2txt[$i]}{$a3[$i][$j]} and do { push(@toadd, $geneID2alnID{$tax2txt[$i]}{$a3[$i][$j]}); next; };
      $is_allAln = 0;
      last;
    }
    $is_allAln == 1 or next;
    for (my $i=0; $i<@toadd; $i++) {
      push(@{$tax2gene[$i]}, $toadd[$i]);
    }
  }
} else {
  @tax2txt = @tax;
  for (my $i=0; $i<@tax2txt; $i++) {
    $tax2gene[$i] = [ map { "$_->[0].$_->[1]" } sort { $a->[1] <=> $b->[1] } @{$alnID_inTax{$tax2txt[$i]}} ];
  }
}

for (my $i=0; $i<@tax2txt; $i++) {
  print STDOUT ">$tax2txt[$i]\n";
  my $ss = "";
  # my $v=-1;
  for my $k1 (@{$tax2gene[$i]}) {
    # $v++;
    unless (defined $s1{$k1}) {
      # warn "i=$i; $tax2txt[$i];\n@{$tax2gene[$i]}[0..5]\n";
      # die "[Err] bad $v-th seq ID [$k1]\n";
      die "[Err] Bad alnID [$k1]\n";
    }
    $ss .= $s1{$k1}{'seq'};
  }
  $ss =~ s!(.{100})!$1\n!g; chomp($ss);
  print STDOUT "$ss\n";
}


