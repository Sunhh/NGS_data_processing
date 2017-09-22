#!/usr/bin/perl 
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 est2genome.psl\n"; 

my @lvls = (0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1); 

my %calc; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	(defined $ta[0] and $ta[0] ne "" and $ta[0] =~ m/^(\d+)$/ and $ta[0] > 0) or next; 
	my $cov_len = $ta[12]-$ta[11]; 
	for my $r1 (@lvls) {
		$cov_len >= $ta[10] * $r1 and $calc{$r1}{$ta[9]}++; 
	}
}
print STDOUT join("\t", qw/Olap_Ratio EST_Number/)."\n"; 
for my $r1 (@lvls) {
	my $num = scalar( keys %{$calc{$r1}} ); 
	print STDOUT join("\t", $r1, $num)."\n"; 
}
