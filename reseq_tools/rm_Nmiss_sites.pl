#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2
	"indel_asN!", # Rule 1 : Treat all indel as N missing; 
	"hete_asN!",  # Rule 2 : Treat heterozygous as missing; 
	"maxMiss:f",  # >= 0.05 
	"only_count!", # Do not do the filtering if given. 
); 

$opts{'maxMiss'} //= 0.05; 
$opts{'startColN'} //= 2; 

my $maxAllowMiss = $opts{'maxMiss'} * 100; 

sub usage {
	print STDERR <<HH; 

perl $0 in.snp > in.snp.filtered

-help 
-startColN        [$opts{'startColN'}]
-maxMiss          [$opts{'maxMiss'}]  missing rate should be <= maxMiss 
-indel_asN        [Bool] If given, genotype with '*|+' will be replaced with N
-hete_asN         [Bool] If given, genotype with m/^[ATGCN]{2,}$/ will be replaced with N 
-only_count       [Bool] Only count the missing rate instead of filter the table. 

HH
	exit(1); 
}

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

while (<>) {
	s/[^\t\S]+$//; 
	my @ta = split(/\t/, $_); 
	my ($chr, $pos) = @ta[0,1]; 
	if ($chr eq 'chr') {
		if ( $opts{'only_count'} ) {
			print STDOUT join("\t", qw/chr pos NmissRate/)."\n"; 
		} else {
			print STDOUT "$_\n"; 
		}
		next; 
	}
	
	# Counting 
	my $missingCnt = 0; 
	my $totalCnt = 0; 
	for my $tb (@ta[$opts{'startColN'}..$#ta]) {
		$tb = uc($tb); 
		$opts{'indel_asN'} and $tb =~ m/[*+]/ and $tb = 'N';  # R1 
		$opts{'hete_asN'} and $tb =~ m/^[ATGCN*]{2,}$/ and $tb = 'N'; # R2 
		$tb eq 'N' and $missingCnt++; 
		$totalCnt ++; 
	}
	my $missingRate = int($missingCnt/$totalCnt*10000+0.5)/100; 
	
	# Filtering
	if ( $opts{'only_count'} ) {
		print STDOUT join("\t", $chr, $pos, $missingRate)."\n"; 
	} else {
		if ( $missingRate <= $maxAllowMiss) {
			print STDOUT "$_\n"; 
		}
	}
}
