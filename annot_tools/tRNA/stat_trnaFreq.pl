#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 arab.chr.fa.trnascan.o.slct > arab.chr.fa.trnascan.o.slct.stat\n"; 

# [Sunhh@bioinfor01 work]$ head -4 arab.chr.fa.trnascan.o.slct
# trna_1  1       +       306384  306456  Val     TAC     ggtgctgtggtgtagtggttatcacgtttgccttacacgcaaaaggtctccagttcgatcctgggcagcacca
# trna_2  1       +       515494  515566  Phe     GAA     gcggggatagctcagttgggagagcgtcagactgaagatctgaaggtcgcgtgttcgatccacgctcaccgca
# trna_3  1       +       552640  552711  His     GTG     gtggctgtagtttagtggtaagaattccacgttgtggccgtggagacctgggctcgaatcccagcagccaca
# trna_4  1       +       604402  604474  Lys     CTT     gcccgtctagctcagttggtagagcgcaaggctcttaaccttgtggtcgtgggttcgagccccacggtgggcg

my %h; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$h{$ta[5]}{$ta[6]}++; 
}
for my $aa ( sort keys %h ) {
	my $sum_gene = 0; 
	for my $cc ( sort keys %{$h{$aa}} ) {
		$sum_gene += $h{$aa}{$cc}; 
	}
}

