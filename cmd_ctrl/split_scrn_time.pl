#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 scrn.t1 | less -S\n"; 

while (<>) {
	chomp; 
	s!^\[([^\[\]]+)\]\[([^\[\]]+)\]!! or next; 
	my ($t1, $t2) = ($1, $2); 
	s!^\s+!!; 
	print join("\t", $t1, $t2, $_)."\n"; 
}
