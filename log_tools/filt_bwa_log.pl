#!/usr/bin/perl
use strict; 
use warnings; 

my @pat_to_rm = ('^\[bwa_aln_core\]', '^\[infer_isize\]', '^\[bwa_sai2sam_pe_core\]', '\[bwa_paired_sw\]'); 

my @pat_use; 
for (@pat_to_rm) {
	push(@pat_use, qr/$_/s );
}

-t and !@ARGV and die "perl $0 log.bwa\nTo skip patterns : @pat_use\n"; 


while (<>) {
	my $is_skip = 0; 
	for my $qPat (@pat_use) {
		m!$qPat! and do { $is_skip = 1; last; }; 
	}
	$is_skip or print; 
}

