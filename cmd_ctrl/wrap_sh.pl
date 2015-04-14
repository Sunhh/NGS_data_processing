#!/usr/bin/perl
use strict; 
use warnings; 

print STDOUT <<HH;
### Basic functions.
function exe_cmd {
	echo \"[\$(date)][CMD] \$1\"
	eval \$1
	echo \"[\$(date)][Rec] Done.\"
}
function tsmsg {
	echo \"[\$(date)]\$1\"
}

HH

while (<>) {
	if (m/^\s*#/ or m/^\s*$/) {
		print; 
	} else {
		s/[^\t \S]+$//; 
		print "exe_cmd \"$_\"\n"; 
	}
}
