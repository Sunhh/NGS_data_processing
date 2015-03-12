#!/usr/bin/perl
use strict; 
use warnings; 

my $head = <>; 
my %cnt; 
while (<>) {
	m/^\s*$/ and next; 
	$cnt{'total'} ++; 
	if ( m/^([ATGC])\1*$/ ) {
		$cnt{"HomoBase"} ++; 
	} elsif ( m/^N+$/i ) {
		$cnt{"N"} ++; 
	} elsif ( m/^\*+$/ ) {
		$cnt{"homoDel"} ++; 
	} elsif ( m/\+/) {
		$cnt{"Ins"} ++; 
	} elsif ( m/^[ATGC][ATGC]$/ ) {
		$cnt{"DiHete"} ++; 
	} else {
		$cnt{"Other"} ++; 
	}
}

for (qw/total N HomoBase homoDel Ins DiHete Other/) {
	print join("\t", $_, $cnt{$_})."\n"; 
}
