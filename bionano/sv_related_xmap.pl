#!/usr/bin/perl 
use strict; 
use warnings; 
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

	my $cov_qLen = abs($ta[4]-$ta[4]); # I don't want add 1 because it shouldn't matter in the current length variations. 
	my $cov_sLen = abs($ta[6]-$ta[5]); 
	my $ttl_qLen = $ta[10]; 
	my $ttl_sLen = $ta[11]; 
	my $is_bad = 1; 
	$cov_qLen >= $opts{'maxCov'} * $ttl_qLen and $is_bad = 0; 
	$cov_sLen >= $opts{'maxCov'} * $ttl_sLen and $is_bad = 0; 

	$is_bad == 1 and print STDOUT "$_\n"; 	
}

