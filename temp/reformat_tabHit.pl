#!/usr/bin/perl
use strict; 
use warnings; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] =~ s!^(\S+) !! or die "$_\n"; 
	$ta[1] = $1; 
	print STDOUT join("\t", @ta)."\n"; 
}
