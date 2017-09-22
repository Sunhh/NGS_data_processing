#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

print STDOUT join("\t", qw/Prefix Raw_RdPairs Kept_RdPairs/)."\n"; 
my %infor; 
while (<>) {
	m!\[Rec\]! or next; 
	if ( m!\[Rec\] There are (\d+) read pairs in total \[([^\[\]\s]+)\]! ) {
		my %tmp; 
		$tmp{total_rdPN} = $1; 
		$tmp{pref} = $2; 
		if ( defined $infor{total_rdPN} and $infor{total_rdPN} ne '' ) {
			print STDOUT join("\t", @infor{qw/pref total_rdPN kept_rdPN/})."\n"; 
			%infor = (); 
		}
		%infor = %tmp; 
	} elsif ( m!\[Rec\] There are .+ (\d+) \([\d.]+\%\) reads kept in both! ) {
		$infor{kept_rdPN} = $1; 
		$infor{kept_perC} = $2; 
	} else {
	} 
}

if ( defined $infor{total_rdPN} and $infor{total_rdPN} ne '' ) {
	print STDOUT join("\t", @infor{qw/pref total_rdPN kept_rdPN/})."\n"; 
	%infor = (); 
}


