#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 

# Rules: 
#  R1. Treat all indel as N missing; 
#  R2. Treat heterozygous not missing; 

my $maxAllowMiss = 0.05 * 100; 
$maxAllowMiss = 0.20 * 100; 

my $is_outRate = 1; 
$is_outRate = 1; 

while (<>) {
	$. % 100e3 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	s/[^\t\S]+$//; 
	my @ta = split(/\t/, $_); 
	my ($chr, $pos, $refB) = @ta[0,1,2]; 
	if ($chr eq 'chr') {
		if ( $is_outRate == 1 ) {
			print STDOUT join("\t", qw/chr pos NmissRate/)."\n"; 
		} else {
			print STDOUT "$_\n"; 
		}
		next; 
	}
	
	# Counting 
	my $missingCnt = 0; 
	my $totalCnt = 0; 
	for my $tb (@ta[3..$#ta]) {
		$tb = uc($tb); 
		if ( length($tb) == 1 ) {
			my @td = &SNP_tbl::dna_d2b($tb); 
			scalar(@td) == 1 or scalar(@td) == 2 or $tb = 'N'; 
		} elsif ( $tb =~ m!^[ATGC][ATGC]$! ) {
			; 
		} else {
			$tb = 'N'; 
		}
		$tb eq 'N' and $missingCnt++; 
		$totalCnt ++; 
	}
	my $missingRate = int($missingCnt/$totalCnt*10000+0.5)/100; 
	
	# Filtering
	if ( $is_outRate == 1 ) {
		print STDOUT join("\t", $chr, $pos, $missingRate)."\n"; 
	} else {
		if ( $missingRate <= $maxAllowMiss) {
			print STDOUT "$_\n"; 
		}
	}
}

