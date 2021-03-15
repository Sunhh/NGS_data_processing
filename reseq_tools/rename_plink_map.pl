#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 plink_fmt.map plink_IDnamed.map\n"; 


my %h; 
while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	$ta[1] eq "." and $ta[1] = "$ta[0]_$ta[3]"; 
	defined $h{$ta[1]} and die "[Err] Repeat site ID [$ta[1]]\n"; 
	$h{$ta[1]} = 1; 
	print join("\t", @ta)."\n"; 
}
