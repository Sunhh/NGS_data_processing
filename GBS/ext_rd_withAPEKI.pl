#!/usr/bin/perl 
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 barcode input.fastq\n\nI use CAGC/CTGC as the required restricted site.\n"; 

my $barC = shift; 
my $addCut1 = "${barC}CAGC"; 
my $addCut2 = "${barC}CTGC"; 
my $ll = length($barC); 

my $n = 0; 
while (my $id = <>) {
	$n ++; 
	$n % 1e6 == 1 and &tsmsg("[Msg] $n reads treated.\n"); 
	my $seq=<>; 
	<>; 
	my $qual = <>; 
	if ($seq =~ m/^(?:$addCut1|$addCut2)/o) {
		my $sub_seq = substr($seq, 0, 64+$ll); 
		$sub_seq =~ m/N/ and next; 
		print STDOUT "$id$seq+\n$qual"; 
	}
}

