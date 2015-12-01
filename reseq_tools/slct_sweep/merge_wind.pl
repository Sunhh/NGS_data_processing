#!/usr/bin/perl
use strict; 
use warnings; 
use mathSunhh; 
my $mat_obj = mathSunhh->new(); 

my %raw_blks; 
my %ord; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ord{$ta[0]} //= $.; 
	push(@{$raw_blks{$ta[0]}}, [$ta[1], $ta[2]]); 
}
my %merged_blks; 
print STDOUT join("\t", qw/chr start end/)."\n"; 
for my $tk (sort {$ord{$a} <=> $ord{$b}} keys %raw_blks) {
	$merged_blks{$tk} = $mat_obj->mergeLocBlk( $raw_blks{$tk} ); 
	for my $tr1 (@{$merged_blks{$tk}}) {
		print STDOUT join("\t", $tk, $tr1->[0], $tr1->[1])."\n"; 
	}
}

