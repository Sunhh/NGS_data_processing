#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

-t and !@ARGV and die "gzip -cd lmyPM_filtV.vcf.gz | perl $0 > lmyPM_filtV_PASS.vcf\n"; 

while (<>) {
	if (m!^#!) {
		print STDOUT $_; 
		next; 
	}
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[6] =~ m!^PASS$!i or next; 
	print STDOUT "$_\n"; 
}

