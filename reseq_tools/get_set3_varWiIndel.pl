#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts,
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
); 

$opts{'startColN'} //= 2; 

sub usage {
	print STDERR <<HH; 

perl $0 in.snp > out_woIndel.snp

-help
-startColN      [$opts{'startColN'}] 
-noHeader       [Bool]

HH
	exit(1); 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

while (<>) {
	$. % 1e6 == 1 and &tsmsg("[Msg] Reading $. lines.\n"); 
	s/[^\S\t]+$//; 
	if ( $. == 1 and !$opts{'noHeader'} ) {
		print "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	my $has_indel = 0; 
	for my $tb (@ta[ $opts{'startColN'} .. $#ta ]) {
		$tb =~ m/\*|\+/ and do { $has_indel = 1; last; }; 
	}
	$has_indel == 1 and print "$_\n"; 
}


