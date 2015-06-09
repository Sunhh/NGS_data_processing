#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

my $geno_col = 3; 
-t and !@ARGV and die "perl $0 in_snp.tbl > in_snp.tbl.set2_varOnlyHete\nPlease note the geno_col=$geno_col\n"; 


while (<>) {
	$. % 1e6 == 1 and &tsmsg("[Msg] Reading $. lines.\n"); 
	s/[^\S\t]+$//; 
	my @ta = split(/\t/, $_); 
	if ( $ta[0] eq 'chr' ) {
		print "$_\n"; 
		next; 
	}
	my $base = 'N'; 
	my $has_diff = 0; 
	for (my $i=$geno_col; $i<@ta; $i++) {
		$ta[$i] = uc($ta[$i]); 
		$ta[$i] =~ m/^[ATGC]$|\*|\+/ or $ta[$i] = 'N'; 
		$ta[$i] eq 'N' and next; 
		$base eq 'N' and $base = $ta[$i]; 
		$base ne $ta[$i] and do { $has_diff = 1; last; }; 
	}
	$has_diff == 0 and print "$_\n"; 
}
