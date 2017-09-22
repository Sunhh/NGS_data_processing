#!/usr/bin/perl
use strict; 
use warnings; 
-t and !@ARGV and die "perl $0 region_list SNP_tbl.cols\n"; 

# region_list 
# ChromID WindS   WindE   D1_MEAN S2_MEAN S3_MEAN W5_MEAN S4_MEAN
# chr1    1       100000  0.3654  0.2870  0.1130  0.4397  0.3642
# chr1    10001   110000  0.3163  0.2158  0.0719  0.4857  0.2919
my $f1 = shift; 
open F,'<',"$f1" or die; 
my %region_list; 
while (<F>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@{$region_list{$ta[0]}}, [@ta[1,2]]); 
}
close F; 

# Merge windows. 
for my $k1 (keys %region_list) {
	my @back_arr; 
	for my $ar ( sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$region_list{$k1}} ) {
		my ($s, $e) = @$ar; 
		if ( scalar(@back_arr) > 0 ) {
			if ( $s <= $back_arr[-1][1] ) {
				$e > $back_arr[-1][1] and $back_arr[-1][1] = $e; 
			} else {
				push(@back_arr, [$s,$e]); 
			}
		} else {
			@back_arr = [ $s, $e ];  
		}
	}
	$region_list{$k1} = \@back_arr; 
}

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[0] eq 'chr' ) {
		print STDOUT "$_\n"; 
		next; 
	}
	my $is_in = 0; 
	for my $ar (@{$region_list{$ta[0]}}) {
		my ($s, $e) = @$ar; 
		$s <= $ta[1] and $e >= $ta[1] and do { $is_in = 1; last; }; 
	}
	$is_in == 1 and print "$_\n"; 
}



