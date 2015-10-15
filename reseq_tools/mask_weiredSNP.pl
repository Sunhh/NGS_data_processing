#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

my $geno_col = 2; 

!@ARGV and die "perl $0 appleSNP_rmClose.cols.set1_woSet2Set3.miss20_heteNotN > appleSNP_rmClose.cols.set1_woSet2Set3.miss20_heteNotN.homoR\nPlease note that the geno_col is $geno_col\n"; 

my $l = <>; 
print $l; 
# print join("\t",qw/ChromID Pos GenoN HomoRatio HeteRatio/)."\n"; 
while (<>) {
	$. % 1e6 == 1 and &tsmsg("line $.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	# my $tot = $#ta- $geno_col +1; 
	#my $tot = 0; 
	#my $homN = 0; 
	#my $hetN = 0; 
	for my $tb ( @ta[$geno_col .. $#ta] ) {
		$tb =~ m/^[ATGC*N]$|^[ATGC]\+[ATGC]+$|^[ATGC*][ATGC*]$/ or $tb = "N"; 
		# Sometimes there will be sth. like 'AG+AAA' heterozygous site, but I don't want it . 
	}
	print join("\t", @ta)."\n"; 
	#my $ratHom = sprintf("%.4f", $homN/$tot) * 100; 
	#my $ratHet = sprintf("%.4f", $hetN/$tot) * 100; 
	#print join("\t", $ta[0], $ta[1], $tot, $ratHom, $ratHet)."\n"; 
}
