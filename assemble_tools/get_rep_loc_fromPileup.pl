#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

# scaffold10      1       A       1       ^].     ?
# scaffold10      2       A       1       .       ?
# scaffold10      3       T       1       .       ?
# scaffold10      4       A       1       .       ?

my $maxDepth = 115 * 2; 
my $maxGap = 1000; 
my %block; 
my @ids; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[3] > $maxDepth ) {
		if ( defined $block{$ta[0]} ) {
			if ( $block{$ta[0]}[-1][1]+1+$maxGap >= $ta[1] ) {
				$block{$ta[0]}[-1][1] = $ta[1]; 
			} else {
				push(@{$block{$ta[0]}}, [$ta[1], $ta[1]]); 
			}
		} else {
			push(@ids, $ta[0]); 
			push(@{$block{$ta[0]}}, [$ta[1], $ta[1]]); 
		}
	}
}

for (@ids) {
	for my $tr1 ( @{$block{$_}} ) {
		print STDOUT join("\t", $_, $tr1->[0], $tr1->[1])."\n"; 
	}
}
