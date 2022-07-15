#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 r1_maker_wiFa.gff3 > r1_maker_woFa.gff3\n"; 

while (<>) {
	m!^\s*#+FASTA\s*$!i and last; 
	print; 
}

