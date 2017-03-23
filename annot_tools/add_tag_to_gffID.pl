#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 tag in.fmt.gff3\n"; 

my $tag = shift; 

while (<>) {
	chomp; 
	if (m/^(\s*$|#)/) {
		print "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	if ($ta[2] =~ m!^(protein_match)$!i) {
		$ta[8] =~ s!^ID=!ID=$tag!; 
	} elsif ($ta[2] =~ m!^(match_part)$!i) {
		$ta[8] =~ s!(^|\s|;)Parent=!$1Parent=$tag!; 
	} elsif ($ta[2] =~ m!^dispersed_repeat$!i) {
		$ta[8] =~ s!^ID=!ID=$tag!; 
	} else {
		die "$_\n"; 
	}
	print join("\t", @ta)."\n"; 
}

