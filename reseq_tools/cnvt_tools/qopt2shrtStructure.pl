#!/usr/bin/perl
use strict; 
use warnings; 

-t !@ARGV and die "perl $0 usetax_var_K_2.qopt > usetax_var_K_2.qopt_f\n"; 


while (<>) {
	chomp; 
	my @ta = split(/\s+/, $_); 
	print join("    ", $., $., "(0)", 1, ":" , @ta)."\n"; 
}


