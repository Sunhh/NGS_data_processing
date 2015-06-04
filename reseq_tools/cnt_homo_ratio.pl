#!/usr/bin/perl 
use strict; 
use warnings; 

my $l = <>; 
print join("\t",qw/ChromID Pos HomoRatio/)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $tot = $#ta-3+1; 
	my $homN = 0; 
	for my $tb ( @ta[3 .. $#ta] ) {
		$tb =~ m/^[ATGC]$/i and $homN ++; 
	}
	my $rat = sprintf("%.4f", $homN/$tot) * 100; 
	print join("\t", $ta[0], $ta[1], $rat)."\n"; 
}
