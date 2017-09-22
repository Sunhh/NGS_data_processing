#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 tag in.fa\n"; 
my $tag = shift; 

while (<>) {
	if (m/^>(\S+)/) {
		$_ = ">$tag.$1 $tag\n"; 
	}
	print; 
}
