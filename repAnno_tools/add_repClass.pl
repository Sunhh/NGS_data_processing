#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 RepClass in.fa\n"; 

my $tag = shift; 
&tsmsg("[Rec] Add [$tag] to repeat name.\n"); 
while (<>) {
	if (m/^\s*>/) {
		s/^>(\S+)/>$1#$tag/ or &stopErr("[Err] $_"); 
	}
	print STDOUT $_; 
}
