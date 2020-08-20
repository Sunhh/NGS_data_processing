#!/usr/bin/perl
use strict; 
use warnings; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[0] eq '#CHROM') {
		print join("\t", @ta[0,1], 'DP_SUM')."\n"; 
		next; 
	}
	my $sumDepth = 0; 
	for my $cnt (@ta[2 .. $#ta]) {
		$cnt eq '.' and next; 
		$sumDepth += $cnt; 
	}
	print join("\t", @ta[0,1], $sumDepth)."\n"; 
}
