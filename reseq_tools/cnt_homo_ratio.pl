#!/usr/bin/perl 
use strict; 
use warnings; 

my $geno_col = 3; 

!@ARGV and die "perl $0 appleSNP_rmClose.cols.set1_woSet2Set3.miss20_heteNotN > appleSNP_rmClose.cols.set1_woSet2Set3.miss20_heteNotN.homoR\nPlease note that the geno_col is $geno_col\n"; 

my $l = <>; 
print join("\t",qw/ChromID Pos HomoRatio/)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $tot = $#ta- $geno_col +1; 
	my $homN = 0; 
	for my $tb ( @ta[$geno_col .. $#ta] ) {
		$tb =~ m/^[ATGC]$/i and $homN ++; 
	}
	my $rat = sprintf("%.4f", $homN/$tot) * 100; 
	print join("\t", $ta[0], $ta[1], $rat)."\n"; 
}
