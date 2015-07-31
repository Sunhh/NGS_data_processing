#!/bin/env perl
use strict; 
use warnings; 

my ($interval_mean, $interval_stdev); 
while (<>) {
	chomp; 
	m!^interval_mean\t([\+\-\d.]+)$! and $interval_mean = $1; 
	m!^interval_stdev\t([\+\-\d.]+)$! and $interval_stdev = $1; 
}
print "INS_mean=$interval_mean\n"; 
print "INS_stdev=$interval_stdev\n"; 
print "INS_cutoff=" . ($interval_mean+3*$interval_stdev) . "\n"; 

