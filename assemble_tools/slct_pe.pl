#!/usr/bin/perl
use strict; 
use warnings; 

my @points_slct = 30e3; 
my $max_allowed = 280e3; 
my $min_allowed = 70e3; 
# For PE points 
@points_slct = (); 
for (my $i=5e3; $i<=70e3; $i+=2e3) {
	push(@points_slct, $i);  
}
# Ref chloroplast
$max_allowed = 80e3; 
$min_allowed = 2e3; 

# For MP points 
#@points_slct = (); 
#for (my $i=30e3; $i<=50e3; $i+=2e3) {
#	push(@points_slct, $i);  
#}
# Ref chloroplast
#$max_allowed = 80e3; 
#$min_allowed = 2e3; 
#@points_slct = (100e3, 130e3, 160e3, 190e3, 220e3); 




while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $is_o = 0; 
	for my $p (@points_slct) {
		if ( $ta[3] <= $p ) {
			$ta[3]+$ta[8] > $p and $ta[3]+$ta[8] < $max_allowed and $ta[3] > $min_allowed and do { $is_o = 1; last; }; 
		} else {
			$ta[3]+$ta[8] < $p and $ta[3] < $max_allowed and $ta[3]+$ta[4] > $min_allowed and do { $is_o = 1; last; }; 
		}
	}
	# $is_o == 1 and print "$_\n"; 
	$is_o == 1 and print "$ta[8]\n"; 
}

