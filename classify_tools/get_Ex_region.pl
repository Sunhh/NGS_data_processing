#!/usr/bin/env perl 
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 P1R02scaf.nt_bn6.MiORGN_join\n"; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($in, $en) = (0,0); 
	if ($ta[6] =~ /In:([\d.]+)/) {
		$in = $1; 
	}
	if ($ta[6] =~ /Ex:([\d.]+)/) {
		$en = $1; 
	}

	$en > 0 or next; 
	
	if ($in == 0) {
		; 
	} elsif ($in <= 3) {
		$en >= $in+3 or next; 
	} else {
		$en >= $in*2 or next; 
	}
	print "$_\n"; 
}
