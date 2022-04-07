#!/usr/bin/perl
# 20180411 Use kill() to monitor PID. 
# [4/7/2022] Invoke 'ps' command instead of kill() function to see if a PID is running, because I need authority to correctly use kill().
use strict; 
use warnings; 
BEGIN {
	use lib "/usr/local/share/perl5/";
	use lib "/Data/Sunhh/src/general/conda/pkgs/perl-parallel-forkmanager-1.17-0/lib/perl5/site_perl/5.22.0/";
}

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
    if ( &if_running($t1) ) {
			$all_done = 0; 
			last; 
		} else {
			delete $pids{$t1}; 
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

sub if_running {
  my ($tp) = @_;
  my @txt = `ps -p $tp -o pid`;
  chomp(@txt);
  my $is_run = 0;
  for my $l1 (@txt) {
    $l1 =~ s!^\+s|\s+$!!g;
    $l1 =~ m!^\s*${tp}\s*$! and do { $is_run = 1; last; };
  }
  return($is_run);
}

