#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

use Getopt::Long;
my %opts;
GetOptions(\%opts, 
	"out:s", 
	"missR_byIndiv:s", 
	"missR_bySite:s", 
	"rm_lmissR:f", # Save values equal to this. 
	"help!"
); 

sub usage {
	my $infor = <<INFOR; 
##################################################
#command:perl $0 <STDIN|parameters>
# 2014-11-12

  -help

# For missing rates 
  -missR_byIndiv   [out_file_name]
  -missR_bySite    [out_file_name]
  -rm_lmissR       [max_miss_rate](0.2) Will save values equal to this. 

##################################################
INFOR
}#End usage()

# Rules: 
#  R1. Treat all indel as N missing; 
#  R2. Treat heterozygous as missing; 

defined $opts{rm_lmissR} or $opts{rm_lmissR} = 0.05; 
my $fh_out = (defined $opts{out}) ? &fileSunhh::openFH($opts{out}, "write") : \*STDOUT; 
defined $opts

my $maxAllowMiss = 0.05 * 100; 
my $maxFiltMiss  = 0.05 * 100; 
$maxAllowMiss = 0.20 * 100; 

my $is_outRate = 1; 
$is_outRate = 1; 

while (<>) {
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
	for my $tb (@ta[3..$#ta]) {
		$tb = uc($tb); 
		$tb =~ m/^[ATGCN]$/ or $tb = 'N'; # R1 & R2
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


