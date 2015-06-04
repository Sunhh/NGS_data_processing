#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 ppid\n"; 

my $ppid = shift; 
my @txt = `ps --ppid $ppid`; 
for (@txt) {
	chomp; 
	m!^\s*(\d+)! or next; 
	my $pid = $1; 
	&exeCmd_1cmd("kill -9 $pid"); 
}

