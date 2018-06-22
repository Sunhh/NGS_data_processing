#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

-t and !@ARGV and die "perl $0 sorted_GVCF_list\n"; 

my $has_header = 0; 

while (my $l = <>) {
	chomp($l); 
	$l =~ m!^\s*($|#)! and next; 
	my ($fn) = (&splitL("\t", $l))[0]; 
	my $ifh = &openFH($fn, '<'); 
	while (<$ifh>) {
		if (m!^#!) {
			$has_header == 1 and next; 
			m!^#CHROM\t! and $has_header = 1; 
			print STDOUT $_; 
			next; 
		}
		print STDOUT $_; 
	}
	close($ifh); 
}
