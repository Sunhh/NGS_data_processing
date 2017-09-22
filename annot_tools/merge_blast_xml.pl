#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
# use mathSunhh; 

-t and !@ARGV and die "perl $0 r9GoodNameSub_00001.fasta.toNr.xml r9GoodNameSub_00002.fasta.toNr.xml\n"; 

my $iter_num = 0; 

my $has_header = 0; 

my @blkType; 
my @blkType_dynam; 

my $is_tail = 0; 
my $txt_tail = ''; 

while (<>) {
	my $is_typeEnd = 0; 
	my $endType; 
	if ( m!^\<([^/]+)\>\s*$! ) {
		push(@blkType, $1); 
		push(@blkType_dynam, $1); 
	} elsif ( m!^\</(.+)\>\s*$! ) {
		$endType = $1; 
		# $endType eq $blkType[-1] or &stopErr("[Err] Bad end of type section. endType=$endType , blkType=$blkType[-1]\n$_\n"); 
		$endType eq $blkType_dynam[-1] or &stopErr("[Err] Bad end of type section dynam. endType=$endType , dynamType=$blkType_dynam[-1]\n"); 
		pop(@blkType_dynam); 
		$is_typeEnd = 1; 
	}

	if (m/^\<\?xml version|^\<\!DOCTYPE BlastOutput/i) {
		# push(@blkType, 'xml'); 
		$has_header or print STDOUT $_; 
		$txt_tail = ''; $is_tail = 0; 
		next; 
	} elsif ( $blkType[-1] eq 'BlastOutput' ) {
		$has_header or print STDOUT $_; 
		$txt_tail = ''; $is_tail = 0; 
		next; 
	} elsif ( $blkType[-1] eq 'BlastOutput_iterations' ) {
		$has_header or print STDOUT $_; 
		$txt_tail = ''; $is_tail = 0; 
		next; 
	} else {
		$has_header = 1; 
	}

	if ( $is_tail == 1 or m!^\</BlastOutput_iterations\>|^\</BlastOutput>! ) {
		$is_tail = 1; 
		$txt_tail .= $_; 
	} elsif ( m!^\s*<Iteration_iter\-num>(\d+)</Iteration_iter\-num>! ) {
		$iter_num ++; 
		s!^(\s*)<Iteration_iter\-num>(\d+)</Iteration_iter\-num>!$1<Iteration_iter\-num>$iter_num</Iteration_iter\-num>!; 
		print STDOUT $_; 
	} elsif (m!^\s*<Iteration_query\-ID>Query_\d+</Iteration_query\-ID>!) {
		s!^(\s*)<Iteration_query\-ID>Query_\d+</Iteration_query\-ID>!$1<Iteration_query\-ID>Query_$iter_num</Iteration_query\-ID>!; 
		print STDOUT $_; 
	} else {
		print STDOUT $_; 
	} 

}

print STDOUT "$txt_tail"; 
