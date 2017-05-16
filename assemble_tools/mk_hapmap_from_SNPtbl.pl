#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
# chr     pos     base    Ref_97103       GS2_RZ900_Comb  GS3_RZ901_Comb  GS4_Sugarlee_Time1      GS5_JX2_Time1   GS6_JLM_Time1   GS7_JXF_Time1   GS8_XHBFGM_Time1
# WM97_Chr01      10      T       T       T       T       N       T       T       T       N       T       T       T       N       N       N       *       N       T
# WM97_Chr01      11      T       T       T       T       N       T       T       T       N       T       T       T       T       N       N       T       N       T

while (<>) {
	chomp; s/[^\S\t]+$//; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		# print STDOUT "$_\n"; 
		print STDOUT join("\t", "rs#", "alleles", 'chrom', "pos", @ta[3..$#ta])."\n"; 
		next; 
	}
	my @dd; 
	my %cnt; 
	for (my $i=3; $i<@ta; $i++) {
		$ta[$i] = uc($ta[$i]); 
		if ( $ta[$i] =~ m/^([ATGC])$/ ) {
			$cnt{allele}{$1} += 2; 
			$ta[$i] .= $ta[$i]; 
		} elsif ( $ta[$i] =~ m/^([ATGC])([ATGC])$/ ) {
			$cnt{allele}{$1} ++; 
			$cnt{allele}{$2} ++; 
		} else {
			$ta[$i] = "NN"; 
			$cnt{N} ++; 
		}
		$cnt{total} ++; 
	}
	for my $t (qw/N/) {
		defined $cnt{$t} or $cnt{$t} = 0; 
	}
	my @alleles = sort { $cnt{allele}{$b} <=> $cnt{allele}{$a} } keys %{ $cnt{allele} }; 
	scalar( @alleles ) == 2 or next; 
	$cnt{N} <= 0.02 * $cnt{total} or next; 
	$ta[0] =~ m/^WM97_Chr0?(\d+)$/ or $ta[0] =~ m/^Chr0?(\d+)$/i or die "$ta[0]\n"; 
	my $chrNum = $1; 
	print STDOUT join("\t", "$ta[0]_$ta[1]", "$alleles[0]/$alleles[1]", $chrNum, $ta[1], @ta[3..$#ta])."\n"; 
}


