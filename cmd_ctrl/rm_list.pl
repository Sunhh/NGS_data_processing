#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

-t and !@ARGV and die "perl $0 filelist_toRM\n"; 

while (<>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	&fileSunhh::_rmtree($ta[0]); 
}

