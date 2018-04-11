#!/usr/bin/perl
# 20180411 Use kill() to monitor PID. 
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"pid:i@", 
	"sleep_seconds:i", # 10
); 

my $help_txt = <<HH; 
################################################################################
#  perl $0 -pid pid_to_wait_1 -pid pid_to_wait_2 ...
#
#  This process will end when none of the '-pid's exists in the system; 
#
#  -help
#
#  -sleep_seconds      [10]
################################################################################
HH

$opts{'sleep_seconds'} //= 10; 

defined $opts{'pid'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
my %pids = map { $_=>1 } @{$opts{'pid'}}; 

while (1) {
	my $all_done = 1; 
	my @tk = sort keys %pids; 
	for my $t1 (@tk) {
		if (kill(0, $t1) == 0) {
			delete $pids{$t1}; 
		} else {
			$all_done = 0; 
			last; 
		}
	}
	$all_done == 1 and last; 
	sleep $opts{'sleep_seconds'}
}


