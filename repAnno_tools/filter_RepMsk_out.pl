#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 in_RepMsk.out\n"; 

while (<>) {
	unless ( m/^\s*\d+/ ) {
		# print; 
		next; 
	}
	chomp; 
	s/^\s+//; s/\s+$//; 
	my @ta = split(/\s+/, $_); 
	my $name1=$ta[4]; 
	my $name2=$ta[9]; 
	$name1 =~ m/^RR\d+_(seq\d+)_(\d+)_(\d+)_INN_(\S+)$/ or die "name1=$name1\n"; 
	my @nn1 = ($1,$2,$3,$4); 
	$name2 =~ m/^([^\s:]+):(\d+)\-(\d+):([FR])$/ or die "$name2\n"; 
	my @nn2 = ($1,$2,$3,$4); 
	my $a = 0; 
	if ( $nn1[3] eq $nn2[0] ) {
		$nn1[1]-1 == $nn2[2] and $a = 1; 
		$nn1[2]+1 == $nn2[1] and $a = 1; 
	}
	$a == 1 or print STDOUT "$_\n"; 
}
