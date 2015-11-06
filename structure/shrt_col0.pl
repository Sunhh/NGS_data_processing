#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 long_indvID_table\n\nI will shorten the first column to 9-10 characters.\n"; 

my $shrt_len = 9; 
my %h; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] = substr($ta[0], 0, $shrt_len); 
	my $tk = $ta[0]; 
	my $suff = "a"; 
	while (defined $h{$tk}) {
		$suff++; 
		$tk = "$ta[0]$suff"; 
	}
	$h{$tk} = 1; 
	$ta[0] = $tk; 
	print join("\t", @ta)."\n"; 
}
