#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 scrn.t1 | less -S\n"; 

while (<>) {
	chomp; 
	# Begin index search ...[Fri Jun 12 10:21:32 2020][CMD_done]gzip clean_files/lnc_GWRT_Rep2_*fq
	if ( s!^\[([^\[\]]+)\]\[([^\[\]]+)\]!! ) {
		my ($t1, $t2) = ($1, $2); 
		s!^\s+!!; 
		print join("\t", $t1, $t2, $_)."\n"; 
	} elsif ( s!^.*\[(\S\S\S +\S\S\S +\d+ +\d+:\d+:\d+ +\d+)\]\[([^\[\]]+)\]!! ) {
		my ($t1, $t2) = ($1, $2); 
		s!^\s+!!; 
		print join("\t", $t1, $t2, $_)."\n"; 
	}
}
