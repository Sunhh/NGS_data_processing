#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 



while (<>) {
	$. % 1e6 == 1 and &tsmsg("[Msg] Reading $. lines.\n"); 
	s/[^\S\t]+$//; 
	if ( m!^chr\t! ) {
		print "$_\n"; 
		next; 
	}
	m/\*|\+/ and print "$_\n"; 
}
