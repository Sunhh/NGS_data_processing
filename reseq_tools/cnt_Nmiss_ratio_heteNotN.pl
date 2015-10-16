#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
); 
$opts{'startColN'} //= 2; 
my $geno_col = $opts{'startColN'}; 

my $help_txt = <<HH; 

Replaced by cnt_homo_hete_ratio.pl script, which can do same thing with some conversion. 

perl $0 in_snp.tbl > in_snp.tbl.missRatio

-help 
-startColN        [$opts{'startColN'}]


Please note geno_col=$geno_col

HH

!@ARGV and -t and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

# Rules: 
#  R1. Treat all indel as N missing; 
#  R2. Treat heterozygous not missing; 


my $maxAllowMiss = 0.05 * 100; 
$maxAllowMiss = 0.20 * 100; 

my $is_outRate = 1; 
$is_outRate = 1; 

while (<>) {
	$. % 100e3 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	s/[^\t\S]+$//; 
	my @ta = split(/\t/, $_); 
	my ($chr, $pos, $refB) = @ta[0,1,2]; 
	if ($chr eq 'chr') {
		if ( $is_outRate == 1 ) {
			print STDOUT join("\t", qw/chr pos NmissRate/)."\n"; 
		} else {
			print STDOUT "$_\n"; 
		}
		next; 
	}
	
	# Counting 
	my $missingCnt = 0; 
	my $totalCnt = 0; 
	for my $tb (@ta[$geno_col .. $#ta]) {
		$tb = uc($tb); 
		if ( length($tb) == 1 ) {
			my @td = &SNP_tbl::dna_d2b($tb); 
			scalar(@td) == 1 or scalar(@td) == 2 or $tb = 'N'; 
		} elsif ( $tb =~ m!^[ATGC][ATGC]$! ) {
			; 
		} else {
			$tb = 'N'; 
		}
		$tb eq 'N' and $missingCnt++; 
		$totalCnt ++; 
	}
	my $missingRate = int($missingCnt/$totalCnt*10000+0.5)/100; 
	
	# Filtering
	if ( $is_outRate == 1 ) {
		print STDOUT join("\t", $chr, $pos, $missingRate)."\n"; 
	} else {
		if ( $missingRate <= $maxAllowMiss) {
			print STDOUT "$_\n"; 
		}
	}
}

