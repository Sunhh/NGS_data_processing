#!/usr/bin/perl
use strict; 
use warnings; 

# ==> ngs.nlis <==
# Key	Length	MatchStart	MatchEnd	MatchLen
# WM97_scaffold11	169061	3288	3336	49
# WM97_scaffold11	169061	14316	14502	187
# WM97_scaffold11	169061	58876	58908	33
# 
# ==> manual_alignments.txt.stat <==
# ID1	WM97pbV0_000000F	16657	113511	96855
# ID1	WM97pbV0_000000F	128093	139235	11143
# ID1	WM97pbV0_000000F	181677	203420	21744
# ID1	WM97pbV0_000000F	207866	208637	772

!@ARGV and die "perl $0 ngs.nlis manual_alignments.txt.stat > ngs.nlis.coverTag\n"; 
my $nlis = shift; 
my $stat = shift; 

open F2,'<',"$stat" or die; 
my %covered; 
while (<F2>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'ID1' and next; 
	push(@{$covered{$ta[1]}}, [@ta[2,3]]); 
}
close F2; 

open F1,'<',"$nlis" or die; 
while (<F1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[0] eq 'Key') {
		print join("\t", @ta, "Covered")."\n"; 
		next; 
	}
	my $is_cover = 0; 
	for my $t1 (@{$covered{$ta[0]}}) {
		$t1->[0] <= $ta[2] and $t1->[1] >= $ta[3] and do { $is_cover = 1; last; };  
	}
	print join("\t", @ta, $is_cover)."\n"; 
}
close F1; 


