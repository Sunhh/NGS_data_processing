#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

my $geno_col = 2; 

sub usage {
	print STDERR <<HH; 

Only simple genotypes /^[ATGC*N]\\$|^[ATGC]\\+[ATGC]+\\$|^[ATGC*][ATGC*]\$/ kept.

perl $0 merged.snp > merged.snp.woWeired

Please note that the geno_col is $geno_col
And first line is not checked. 

HH
	exit(1); 
}

!@ARGV and &usage(); 

my $l = <>; 
print $l; 
# print join("\t",qw/ChromID Pos GenoN HomoRatio HeteRatio/)."\n"; 
while (<>) {
	$. % 1e6 == 1 and &tsmsg("line $.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	for my $tb ( @ta[$geno_col .. $#ta] ) {
		$tb =~ m/^[ATGC*N]$|^[ATGC]\+[ATGC]+$|^[ATGC*][ATGC*]$/ or $tb = "N"; 
		# Sometimes there will be sth. like 'AG+AAA' heterozygous site, but I don't want it . 
	}
	print join("\t", @ta)."\n"; 
}
