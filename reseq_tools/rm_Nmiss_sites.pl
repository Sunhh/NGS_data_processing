#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

# Rules: 
#  R1. Treat all indel as N missing; 
#  R2. Treat heterozygous as missing; 

my $maxAllowMiss = 0.05 * 100; 
my $maxFiltMiss  = 0.05 * 100; 
$maxAllowMiss = 0.05 * 100; 

my $is_outRate = 1; 
$is_outRate = 0; 

while (<>) {
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
		$tb =~ m/^[ATGCN]$/ or $tb = 'N'; # R1 & R2
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
