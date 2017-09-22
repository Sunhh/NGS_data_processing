#!/usr/bin/perl
# 2016-02-18 Further trimming '(Similarity to unknown protein)'; 
use strict; 
use warnings; 
use fileSunhh; 
while (<>) {
	if ( m/^\s*(#|$)/ ) {
		next; 
	}
	chomp; 
	my @ta = &splitL("\t", $_); 
	if ( defined $ta[3] ) {
		my $prev_ss = $ta[3]; 
		my $aftr_ss = &trim($prev_ss); 
		while ( $prev_ss ne $aftr_ss ) {
			$prev_ss = $aftr_ss; 
			$aftr_ss = &trim($prev_ss); 
		}
		$ta[3] = $aftr_ss; 
	}
	print join("\t", @ta)."\n"; 
}

sub trim {
	my $sss = shift; 
	$sss =~ s!\(\s+!(!g; $sss =~ s!\s+\)!)!g; 
	$sss =~ s!^([^()]*?)\s*\)!$1!; 
	$sss =~ s!^\s+|\s+$!!g; 
	$sss =~ s!\{\s*(ECO:\d+\|[^\s:]+:[^\s:]+(,\s*)*?)*\s*\}!!g; 
	$sss =~ s!\{\s*ECO:\d+\s*\}!!g; 
	$sss =~ s!\(protein\)!!ig; 
	$sss =~ s!\(putative [^()]+\)!!ig; $sss =~ s!^\s+|\s+$!!g; 
	$sss =~ s!\(Uncharacterized protein\)!!ig; $sss =~ s!^\s+|\s+$!!; 
	$sss =~ s!^\(([^()]+)\)$!$1!g; 
	$sss =~ s!^\(([^()]*\([^()]*\)[^()]*)\)$!$1!g; 
	$sss =~ s!\(\s*\)!!g; 
	$sss =~ s!^([^()]+?)\s*\)\s*$!$1!; 
	$sss =~ s!^([^()]*\([^()]+\)[^()]*)\s*\)\s*$!$1!; 
	$sss =~ s!^\s+|\s+$!!; 
	$sss =~ s!^protein$!!i; 
	$sss =~ s!\(Similarity to unknown protein\)!!i; 
	$sss =~ s!^Full=!!; 

	$sss =~ s!^\s+$!!; 
	$sss eq '' and $sss = "Unknown protein"; 
	return $sss; 
}
