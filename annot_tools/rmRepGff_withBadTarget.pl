#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

my $help = <<HH; 
perl $0 ID_class.not4Annot wm97pb_v2ID.scf.fa.mskAll_4maker.gff > wm97pb_v2ID.scf.fa.mskAll_4maker_clean.gff

# head -4 ID_class.not4Annot
rnd-1_family-53 Simple_repeat
rnd-1_family-410        Simple_repeat
rnd-5_family-2673       Simple_repeat
rnd-5_family-2933       Simple_repeat

HH

!@ARGV and die $help; 

my $f1 = shift; 
my %badID = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $f1 ); 

# Example of repeat_4maker.gff 
# ##gff-version 3
# ClaScf_0001     RepeatMasker    dispersed_repeat        3       6175    7348    -       .       ID=0:1;Target=RR335_seq19_2552196_2555970_LTR_ClaScf_0020 1495 2858
# ClaScf_0001     RepeatMasker    dispersed_repeat        6173    8771    22262   -       .       ID=0:2;Target=RR100_seq1_4243177_4254363_LTR_ClaScf_0002 4360 6957
# ClaScf_0001     RepeatMasker    dispersed_repeat        8747    8928    883     -       .       ID=0:3;Target=RR794_seq6_15451130_15461819_LTR_ClaScf_0007 3392 3578
while (<>) {
	m!^\s*(#|$)! and do { print; next; }; 
	chomp; 
	my @ta=split(/\t/, $_); 
	$ta[8] =~ m!;Target=(\S+)! or die "Bad ta[8]: $_\n"; 
	defined $badID{$1} and next; 
	print "$_\n"; 
}


