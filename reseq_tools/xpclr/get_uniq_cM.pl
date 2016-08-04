#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 

!@ARGV and die "perl $0   marker_loc2GM_man\n"; 

my $fn_gm  = shift; 

my $fh_gm  = &openFH($fn_gm, '<'); 
my %gmP_to_line; 
my %ord; 
while (&wantLineC($fh_gm)) {
	my @ta = &splitL("\t", $_); 
	if ( $. == 0 and $ta[0] eq 'MarkerID' ) {
		print STDOUT "$_\n"; 
		next; 
	}
	my $tk = "$ta[7]_$ta[8]"; 
	push(@{$gmP_to_line{$tk}}, [@ta]); 
	$ord{$tk} //= $.; 
}
close($fh_gm); 

for my $tk (sort { $ord{$a} <=> $ord{$b} } keys %gmP_to_line) {
	my $midN = int( $#{$gmP_to_line{$tk}}/2 ); 
	print STDOUT join("\t", @{$gmP_to_line{$tk}[$midN]})."\n"; 
}

