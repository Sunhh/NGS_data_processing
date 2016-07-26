#!/usr/bin/perl
use strict; 
use warnings; 

my %prev; 
my @lines; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'MarkerID' and next; 
	$ta[1] eq $ta[7] or die "Bad:$_\n"; 
	push(@lines, [$_, '']); 
	my %curr; 
	$curr{'chrID'} = $ta[1]; 
	$curr{'chrP'}  = $ta[2]; 
	$curr{'chrMend'} = $ta[8]; 
	$curr{'chrLine'} = $_; 

	if ( !(defined $prev{'chrID'}) or $prev{'chrID'} ne $curr{'chrID'} ) {
		%prev = %curr; 
		next; 
	}
	
	if ( $prev{'chrMend'} < $curr{'chrMend'} and $prev{'chrP'} > $curr{'chrP'} ) {
		$lines[-2][1] = 'Chk:'; 
		$lines[-1][1] = 'Chk:'; 
	}

	%prev = %curr; 
}
for (@lines) {
	print "$_->[1]$_->[0]\n"; 
}
