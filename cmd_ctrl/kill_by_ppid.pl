#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 ppid curr_user\n"; 

my $ppid = shift; 
my $curr = shift; 
$curr //= 0; 

my @txt = `ps --ppid $ppid`; 
my %cur; 
if ( $curr ) {
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

