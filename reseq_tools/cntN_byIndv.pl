#!/usr/bin/perl
use strict; 
use warnings; 


use LogInforSunhh; 
my $head = <>; 
chomp($head); 
my @ha=split(/\t/, $head); 
my @cnt; 
while (<>) { 
	($.-1) % 1e4 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	for (my $i=2; $i<@ta; $i++) { 
		$ta[$i]=~m/^[ATGC]$/ and $cnt[$i]++; 
	} 
} 
for (my $i=0; $i<@ha; $i++) { 
	$cnt[$i]//=0; 
	print "$ha[$i]\t$cnt[$i]\n";
}


