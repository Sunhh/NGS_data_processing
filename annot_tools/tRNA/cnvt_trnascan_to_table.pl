#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

my %h; 
while (<>) {
	chomp; 
	if (m!^sequence name=\s*(\S+)!) {
		$h{'seqID'} = $1; 
	} elsif (m!^start position!) {
		m!^start position=\s*(\d+) end position=\s*(\d+)! or die "$_\n"; 
		$h{'eleN'} ++; 
		$h{'eleS'} = $1; 
		$h{'eleE'} = $2; 
		$h{'eleStr'} = '+'; 
		if ( $h{'eleS'} > $h{'eleE'} ) {
			$h{'eleStr'} = '-'; 
			($h{'eleS'}, $h{'eleE'}) = ($h{'eleE'}, $h{'eleS'}); 
		}
	} elsif (m!^potential tRNA sequence=\s*(\S+)!) {
		$h{'eleSeq'} = $1; 
	} elsif (m!^(D signal=|amino-acyl stem=|D stem=|anticodon stem=|TpsyC stem=)!) {
		; 
	} elsif (m!^tRNA predict as a tRNA!) {
		m!^tRNA predict as a tRNA\-\s*(\S{3}) : anticodon (\S{3})! or die "$_\n"; 
		$h{'eleType'} = $1; 
		$h{'eleAntiCodon'} = uc($2); 
		print STDOUT join("\t", 
			"trna_$h{'eleN'}", 
			@h{qw/seqID eleStr eleS eleE eleType eleAntiCodon eleSeq/}
		)."\n"; 
		for my $t0 (qw/eleStr eleS eleE eleType eleAntiCodon eleSeq/) {
			$h{$t0} = ''; 
		}
	} elsif (m!^number of base pairing in the anticodon stem|^\s*$!) {
		for my $t0 (qw/eleStr eleS eleE eleType eleAntiCodon eleSeq/) {
			$h{$t0} = ''; 
		}
	} elsif (m!^potential intron between positions !) {
		; 
	} elsif (m!^complementary strand!) {
		; 
	} elsif (m!^anticodon includes unknown bases!) {
		; 
	} elsif (m!^number of !) {
		; 
	} else {
		die "unknown line: $_\n"; 
	}
}


