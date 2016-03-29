#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts,
	"help!", 
	"forAlluser!", 
	"doWith_1!", 
); 

my $help_txt = <<HH; 

perl $0 ppid_to_Kill 

-help

-forAlluser       [Boolean] If given, try to kill all child processes of input ppid. If not given, only consider within current session. 

-doWith_1         [Boolean] Have to be given if ppid is 1. This is dangerous. 

HH

my $help_1 = <<AA;

Input ppid is 1, please add -doWith_1 for continue. 

AA

!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
my $ppid = shift; 
$ppid == 1 and &LogInforSunhh::usage($help_1); 


my @txt = `ps --ppid $ppid`; 

my %cur; 
unless ( $opts{'forAlluser'} ) {
	my @txt_cur = `ps`; 
	for (@txt_cur) {
		chomp; 
		m!^\s*(\d+)! or next; 
		my $pid = $1; 
		$cur{'cur_pid'}{$pid} = 1; 
	}
}

for (@txt) {
	chomp; 
	m!^\s*(\d+)! or next; 
	my $pid = $1; 
	defined $cur{'cur_pid'}{$pid} or next; 
	&exeCmd_1cmd("kill -9 $pid"); 
}

