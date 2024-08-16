#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "cat 02_ma_mo_byScf/ma_mo_byChr.coll.ks.tab.mo_mo | perl $0 > data/dup.txt\n"; 

while (<>) { 
	chomp; 
	my @ta=split(/\t/, $_); 
	@ta = @ta[1,2,3,4,5,6,7]; 
	$ta[0] eq 'Chrom1' and next; 
	$ta[1]>$ta[2] and @ta[1,2]=@ta[2,1]; 
	$ta[4]>$ta[5] and @ta[4,5]=@ta[5,4]; 
	$ta[6] eq "-" and @ta[5,4]=@ta[4,5]; 
	$ta[1]--; 
	$ta[4]--; 
	print join("\t", @ta[0..5] )."\n"; 
}

