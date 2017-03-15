#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"gapLen:i", # 1000
); 
$opts{'gapLen'} //= 1000; 

-t and !@ARGV and die "perl $0 -gapLen 1000 in.agp > o.agp\n"; 

my %aa; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$aa{$ta[0]}{'ord'} //= $.; 
	push(@{$aa{$ta[0]}{'arr'}}, [@ta]); 
}
for my $k1 (sort { $aa{$a}{'ord'} <=> $aa{$b}{'ord'} } keys %aa) {
	my @ta2; 
	for my $t1 (@{$aa{$k1}{'arr'}}) {
		print join("\t", @$t1)."\n"; 
		@ta2 = @$t1; 
	}
	print join("\t", $ta2[0], $ta2[2]+1, $ta2[2]+$opts{'gapLen'}, $ta2[3]+1, 'N', $opts{'gapLen'}, 'scaffold', 'yes', 'map')."\n";
	print join("\t", $ta2[0], $ta2[2]+$opts{'gapLen'}+1, $ta2[2]+$opts{'gapLen'}+1, $ta2[3]+2, 'W', 'NA', 1,1, '+')."\n";  
}
