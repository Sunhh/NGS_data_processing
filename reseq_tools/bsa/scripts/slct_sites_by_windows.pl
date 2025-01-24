#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

!@ARGV and die "perl $0 position_list  site.al\n"; 

my $fn1 = shift; # windows; 
my $fn2 = shift; # sties; 

my $fh1 = &openFH($fn1, '<'); 
my %required; 
while (<$fh1>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	push(@{$required{$ta[0]}}, [@ta[1,2]]); 
}
close($fh1); 
for my $k1 (keys %required) {
	@{$required{$k1}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$required{$k1}}; 
}

my $fh2 = &openFH($fn2, '<'); 
while (<$fh2>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	if ($. == 1 and $ta[0] =~ m!^(chrom$|chr$)!i and $ta[1] =~ m!^(pos$|position$)!i) {
		print STDOUT "$_\n"; 
		next; 
	}
	defined $required{$ta[0]} or next; 
	my $is=0; 
	for my $a1 (@{$required{$ta[0]}}) {
		$a1->[0] > $ta[1] and last; 
		$a1->[1] < $ta[1] and next; 
		$is = 1; 
		last; 
	}
	$is == 1 and print STDOUT "$_\n"; 
	
}
close($fh2); 

