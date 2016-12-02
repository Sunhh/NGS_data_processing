#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

-t and !@ARGV and &LogInforSunhh::usage("\nperl $0 sample.1col > sample.1col.typeC\n\n"); 

my $fh = \*STDIN; 

@ARGV > 0 and $fh = &openFH($ARGV[0], '<'); 

my $head = <$fh>; 
my %cnt; 
while (<$fh>) {
	m/^\s*$/ and next; 
	$cnt{'total'} ++; 
	if ( m/^([ATGC])\1*$/ ) {
		$cnt{"HomoBase"} ++; 
	} elsif ( m/^N+$/i ) {
		$cnt{"N"} ++; 
	} elsif ( m/^[ATGC][ATGC]$/ ) {
		$cnt{"DiHete"} ++; 
	} elsif ( m/^[ATGC\*][ATGC\*]$/ ) {
		$cnt{'heteDel'} ++; 
	} elsif ( m/^\*+$/ ) {
		$cnt{"homoDel"} ++; 
	} elsif ( m/^[ATGCN]\+/) {
		$cnt{"homoIns"} ++; 
	} elsif ( m/^[ATGCN*][ATGCN]+/ ) {
		$cnt{'heteIns'} ++; 
	} else {
		$cnt{"Other"} ++; 
	}
}

chomp($head); 
print join("\t", 'Type', $head)."\n"; 
for (qw/total N HomoBase homoDel heteDel homoIns heteIns DiHete Other/) {
	print join("\t", $_, $cnt{$_})."\n"; 
}
