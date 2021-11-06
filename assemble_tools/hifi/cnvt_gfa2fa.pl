#!/usr/bin/perl
# Convert hifiasm.p_ctg.gfa file to fasta file. 
use strict; 
use warnings; 
# awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa  # get primary contigs in FASTA
!@ARGV and -t and die "perl $0 hifiasm_asm.p_ctg.gfa > hifiasm_asm.p_ctg.fa\n# Used to get the primary contigs in FASTA format.\n\n"; 

while (<>) {
	chomp; 
	m!^S\s+(\S+)\s+(\S+)! or next; 
	print ">$1\n$2\n"; 
}

