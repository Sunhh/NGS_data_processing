#!/usr/bin/perl 
use strict; 
use warnings; 
use mathSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"minConfidence:f", # 0 
	"maxCov:f", # 0 [0-1]
); 

$opts{'minConfidence'} //= 0; 
$opts{'maxCov'} //= 0; 

my $help_txt = <<HH; 

perl $0  -minConfidence    $opts{'minConfidence'}       -maxCov     $opts{'maxCov'}     query_to_anchor.xmap 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

$opts{'maxCov'} > 1 and do { &tsmsg( "[Wrn] -maxCov should not be bigger than 1.\n" ); $opts{'maxCov'} = 1; }; 



while (<>) {
	if (s/^#h XmapEntryID\t/#XmapEntryID\t/) {
		print STDOUT "$_"; 
	}
	m/^\s*(#|$)/ and next; 
	chomp; 
	my @ta = split(/\t/, $_); 

	$ta[8] >= $opts{'minConfidence'} or next; 

	my $str      = $ta[7]; 
	my $cov_qLen = abs($ta[4]-$ta[3]); # I don't want add 1 because it shouldn't matter in the current length variations. 
	my $cov_sLen = abs($ta[6]-$ta[5]); 
	my $ttl_qLen = $ta[10]; 
	my $ttl_sLen = $ta[11]; 
	my $max_cov_qLen = ( $str eq '-' ) ? 
	    &mathSunhh::min($ttl_qLen-$ta[3], $ta[5]-1) + $cov_qLen + &mathSunhh::min($ta[4]-1, $ttl_sLen-$ta[6]-1)  
	  : &mathSunhh::min($ta[3]-1, $ta[5]-1) + $cov_qLen + &mathSunhh::min( $ttl_qLen-$ta[4]-1, $ttl_sLen-$ta[6]-1 ) ; 
	my $max_cov_sLen = ( $str eq '-' ) ? 
	    &mathSunhh::min($ttl_qLen-$ta[3], $ta[5]-1) + $cov_sLen + &mathSunhh::min($ta[4]-1, $ttl_sLen-$ta[6]-1)  
	  : &mathSunhh::min($ta[3]-1, $ta[5]-1) + $cov_sLen + &mathSunhh::min( $ttl_qLen-$ta[4]-1, $ttl_sLen-$ta[6]-1 ) ; 
	# my $max_cov_qLen = &mathSunhh::min( &mathSunhh::max(@ta[3,4]), $ttl_qLen-&mathSunhh::min(@ta[3,4]) ); 
	# my $max_cov_sLen = &mathSunhh::min( &mathSunhh::max(@ta[5,6]), $ttl_qLen-&mathSunhh::min(@ta[5,6]) ); 
	my $is_bad = 1; 
	$cov_qLen >= $opts{'maxCov'} * $max_cov_qLen and $is_bad = 0; 
	$cov_sLen >= $opts{'maxCov'} * $max_cov_sLen and $is_bad = 0; 

	$is_bad == 1 and print STDOUT "$_\n"; 	
}

