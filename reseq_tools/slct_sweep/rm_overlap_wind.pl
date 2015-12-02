#!/usr/bin/perl
use strict; 
use warnings; 

# ChromID WindS   WindE   WindL   BpCnt   perKb_0
# SpoScf_00001    1       10000   10000   10000   0.0924198651487334
# SpoScf_00001    1001    11000   10000   10000   0.0439350166638849
# SpoScf_00001    2001    12000   10000   10000   0.0348612437039331

my @prev; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if (!@prev or $prev[0] ne $ta[0]) {
		print "$_\n"; 
		@prev = @ta; 
		next; 
	}
	$prev[2] >= $ta[1] and next; 
	print "$_\n"; 
	@prev = @ta; 
}

