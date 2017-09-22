#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

# Rules: 
#  R1. Remove adjacent SNP sites within 5bp. 

my $within_dist = 5; 

my %prev; 
while (<>) {
	s/[^\t\S]+$//; 
	my @ta = split(/\t/, $_); 
	my ($chr, $pos, $refB) = @ta[0,1,2]; 
	if ($chr eq 'chr') {
		print STDOUT "$_\n"; 
		next; 
	}

	# Rule 1: 
	my %curr; 
	$curr{'chr'} = $chr; 
	$curr{'pos'} = $pos; 
	$curr{line} = $_; 
	$curr{is_good} = 1; 
	if (scalar(keys %prev) == 0 or $prev{chr} ne $chr) {
		defined $prev{'is_good'} and $prev{'is_good'} == 1 and print STDOUT "$prev{'line'}\n"; 
	} else {
		my $dist2prev = $pos - $prev{pos}+1; 
		if ( $dist2prev <= $within_dist ) {
			# Both are bad. 
			$curr{is_good} = 0; 
			# $prev{is_good} = 0; 
		} else {
			$prev{is_good} == 1 and print STDOUT "$prev{line}\n"; 
		}
	}
	%prev = %curr; 
}
$prev{is_good} == 1 and print STDOUT "$prev{line}\n"; 
