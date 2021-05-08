#!/usr/bin/perl
use strict; 
use warnings; 

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	if ($. == 1) {
		print join("\t", @ta, "absKME")."\n"; 
	} else {
		print join("\t", @ta, abs($ta[1]))."\n"; 
	}
}
