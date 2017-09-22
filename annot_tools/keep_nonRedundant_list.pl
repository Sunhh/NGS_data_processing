#!/usr/bin/perl
use strict;
use warnings; 

my %h; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $tk = "$ta[0]\t$ta[1]"; 
	$h{$tk} = [ @ta[12,13,10] ]; 
}

my %skip; 
my %rm; 
for my $tk ( sort {$h{$a}[2] <=> $h{$b}[2] || $h{$b}[0] <=> $h{$a}[0]} keys %h ) {
	my ($id1, $id2) = split(/\t/, $tk); 
	my ($len1, $len2, $score) = @{$h{$tk}}; 
	defined $skip{$tk} and next; 
	if ( $len2 > $len1 ) {
		if ( !defined $rm{$id2} ) {
			$rm{$id1} = $h{$tk}[0]; 
		}
	} else {
		if ( !defined $rm{$id1} ) {
			$rm{$id2} = $h{$tk}[1]; 
		}
	}
	$skip{"$tk"} = 1; 
	$skip{"$id2\t$id1"} = 1; 
}
for (sort keys %rm) {
	print "$_\t$rm{$_}\n"; 
}

