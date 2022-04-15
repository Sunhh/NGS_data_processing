#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;

-t and !@ARGV and die "perl $0 pref_list\n";

my $pl_trim2CDS_gff = 'perl tools/fmt_gff_trim2CDS.pl';
my $fn_chrList      = 'input/map.CL_CM_CA_CC_WCG_Kord.tab';
my $pl_dealTab      = 'deal_table.pl';
my $pl_dealFas      = 'deal_fasta.pl';

-e "input" or mkdir("input/");
-e "input/map_chr_list" or mkdir("input/map_chr_list/");

my @prefs;
while (<>) {
  chomp;
  # Shrink gene features to CDS-only region.
  &runCmd("$pl_trim2CDS_gff common_data/updated_annotation/${_}.gff3 > input/${_}.trim2CDS.gff3");
  push(@prefs, $_);
}

# Make chr mapping list.
my %pref2col = qw(
ASM1002_97103      0
ASM1004_USVL531    1
ASM1005_PI482246   2
ASM1003_PI537277   3
WCGv2              4
cordophanusV1.5    5
);
for (my $i1=0; $i1<@prefs; $i1++) {
  for (my $i2=$i1+1; $i2<@prefs; $i2++) {
    open F1,'<',"$fn_chrList" or die;
    open O1,'>',"input/map_chr_list/map.$prefs[$i1].vs.$prefs[$i2]" or die;
    open O2,'>',"input/map_chr_list/map.$prefs[$i2].vs.$prefs[$i1]" or die;
    while (<F1>) {
      chomp;
      my @ta=split(/\t/, $_);
      print O1 join(",", @ta[ $pref2col{ $prefs[$i1] }, $pref2col{ $prefs[$i2] } ])."\n";
      print O2 join(",", @ta[ $pref2col{ $prefs[$i2] }, $pref2col{ $prefs[$i1] } ])."\n";
    }
    close O1;
    close O2;
    close F1;
  }
  &runCmd("$pl_dealFas -attr key common_data/$prefs[$i1].ref.fa | tail -n +2 | $pl_dealTab -kSrch_drop -kSrch_idx $fn_chrList -kSrch_idxCol $pref2col{$prefs[$i1]} > input/map_chr_list/unplaced_seqID.$prefs[$i1]");
}

