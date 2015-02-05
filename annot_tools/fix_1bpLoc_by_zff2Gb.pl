#!/usr/bin/perl
use strict; 
use warnings; 

while (<>) {
	while (s/complement\((\d+)\)/$1..$1/) {
		1; 
	}
	print; 
}
