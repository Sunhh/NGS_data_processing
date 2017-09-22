#!/usr/bin/perl
# 2015-06-29 Extend usage of $min2 to col1. 
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 col_grp1 col_grp2 min2_denom in_PIperKb.avg\n"; 

my $col1 = shift; 
my $col2 = shift; 
my $min2 = shift; 
$min2 //= 0; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $v = -1; 
	$. == 1 and $ta[0] =~ m/^(ChrID|ChromID)$/i and do { $v = "C${col1}_to_C${col2}_min${min2}"; print "$v\n"; next; }; 
	($ta[$col2] < $min2 and $ta[$col2] >= 0) and $ta[$col2] = $min2; 
	($ta[$col1] < $min2 and $ta[$col1] >= 0) and $ta[$col1] = $min2; # Added 2015-06-29 
	if ( $ta[$col2] >= $min2 and $ta[$col2] > 0 ) {
		$v = $ta[$col1] / $ta[$col2] ; 
	}
	print "$v\n"; 
}

