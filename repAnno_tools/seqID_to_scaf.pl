#!/usr/bin/perl
use strict;
use warnings; 

my %h; 
print STDERR join("\t", qw/seqID scfID/)."\n"; 
while (<>) {
	if ( m/^>/ ) {
		m/^>(\S+) \(dbseq\-nr (\d+)\) \[(\d+),(\d+)\]$/ or die "$_\n"; 
		my $seqID = "seq$2"; 
		my $scfID = "$1"; 
		my ($eleS, $eleE) = ($3,$4); 
		$_ = ">${seqID}_${eleS}_${eleE}_$scfID\n"; 
		if ( defined $h{$seqID} ) {
			$h{$seqID} eq $scfID or die "$h{$seqID} eq $scfID\n$_\n"; 
		} else {
			print STDERR join("\t", $seqID, $scfID)."\n"; 
			$h{$seqID} = $scfID; 
		}
	}
	print; 
}
