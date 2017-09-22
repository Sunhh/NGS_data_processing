#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"pid:i", 
	"ps_o:s", # ppid,stime,etime,command
	"sleep_seconds:i", # 10
); 

my $help_txt = <<HH; 
################################################################################
#  perl $0 -pid pid_to_wait 
#
#  -help
#
#  -ps_o               ['ppid,stime,etime,command']
#  -sleep_seconds      [10]
################################################################################
HH

$opts{'sleep_seconds'} //= 10; 
$opts{'ps_o'} //= 'ppid,stime,etime,command'; 

defined $opts{'pid'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

while (1) {
	my @aa = `ps -p $opts{'pid'} -o pid,$opts{'ps_o'}`; 
	if (defined $aa[1]) {
		chomp($aa[1]); 
		$aa[1] =~ s!^\s+|\s+$!!g; 
		my @bb = split(/\s+/, $aa[1]); 
		$bb[0] == $opts{'pid'} or &stopErr("[Err] Failed to get PID [$opts{'pid'}]\n"); 
		&tsmsg("[Msg] Waiting for process : $aa[1]\n"); 
		sleep $opts{'sleep_seconds'}; 
		next; 
	}
	last; 
}


