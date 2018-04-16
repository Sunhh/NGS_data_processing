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
	"cmd:s@", 
); 

my $help_txt = <<HH; 
################################################################################
#  perl $0 -pid pid_to_wait_1 -pid pid_to_wait_2 ...
#
#  This process will end when none of the '-pid's exists in the system; 
#
#  -help
#
#  -cmd    [string\@] Excute commands if given. 
#
#  -sleep_seconds      [10]
################################################################################
HH

$opts{'sleep_seconds'} //= 10; 

defined $opts{'pid'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
my %pids = map { $_=>1 } @{$opts{'pid'}}; 
my @cmd; 
defined $opts{'cmd'} and @cmd = @{$opts{'cmd'}}; 

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
	if ($all_done == 1) {
		for my $c (@cmd) {
			&exeCmd_1cmd($c); 
		}
		last; 
	}
	sleep $opts{'sleep_seconds'}
}


