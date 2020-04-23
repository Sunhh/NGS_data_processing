#!/usr/bin/perl
use strict; 
use warnings; 

my %loc; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@{$loc{$ta[1]}}, [@ta[0,2,3]]); 

}


for my $k ( sort keys %loc ) {

	my @loc = sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } @{$loc{$k}}; 	
	for (my $i=0; $i<@loc; $i++){
		my $has = 0; 
		for (my $j=$i+1; $j<@loc; $j++) {
			if ($loc[$i][2] < $loc[$j][1] ) {
				last;  
			} elsif ( $loc[$i][2] >= $loc[$j][1] and $loc[$i][1] <= $loc[$j][2] ) {
				print "$loc[$i][0]\t$loc[$j][0]\n"; 
				$has = 1; 
			} else {
				die "@{$loc[$i]}\n@{$loc[$j]}\n"; 
			}
		}
		if ($has == 0) {
			print "$loc[$i][0]\t\n"; 
		}
	}
} 



