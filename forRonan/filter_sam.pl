#!/usr/bin/perl -w
# 2015-07-14 For ronan's internship project. 
#  Filter parent_Rd to parent_Asm alignments
#   1. Only unique alignments kept. 
#   2. Only 100% match alignments kept. 
#   3. Only mapping quality >= 1 kept. 
#  After filtering, only alignments exactly the same as the reference are kept, which can be used to call same region between Rd_reads and Asm_reference. 
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 

-t and !@ARGV and die "perl $0 in_rd2Asm.sam\n"; 

while (<>) {
	m!^\@! and do { print; next; }; 
	m!\t(?:XT:A:U|NM:i:0)(?:\t|$)! or next; 
	my @ta = split(/\t/, $_); 
	$ta[4] >= 1 or next; 
	print ; 
}



