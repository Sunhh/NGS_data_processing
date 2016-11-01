#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and &LogInforSunhh::usage("perl $0 in_sMao.snp > in_sMao.snp.tab\n"); 

while (<>) { 
	chomp; 
	my @ta=split(/\t/, $_); 
	if ( $. == 1 ) {
		print join("\t", @ta)."\n"; 
		next; 
	}
	for my $tb (@ta[3..$#ta]) { 
		if      ( $tb =~ s!^([ATGCN*])$!$1/$1! ) {
			; 
		} elsif ( $tb =~ s!^([ATGC*])([ATGC*])$!$1/$2! ) {
			; 
		} elsif ( $tb =~ s!^([ATGC])\+([ATGCN]+)$!$1$2/$1$2!) {
			;
		} elsif ( $tb =~ s!^([ATGC])([ATGC])\+([ATGC]+)$!$1/$2$3! ) {
			; 
		} elsif ( $tb =~ m!^\+!) { 
			$tb = "./."; 
		} else { 
			die "$tb\n"; 
		}
		$tb=join("/", sort ($tb =~ m!^([^/]+)/([^/]+)$!) );  
	} 
	print join("\t", @ta)."\n";  
}

