#!/usr/bin/perl
use strict; 
use warnings; 

while (<>) {
	m!^\s*$! and next; 
	m!^\s*#! and next; 
	chomp; 
	my @ta=split(/\t/, $_); 
	if ($. == 1) {
		print join("\t", @ta, qw/mrnaID SeqID MotifStart MotifEnd MotifStr/)."\n"; 
		next; 
	}
	# Parse ID;
	$ta[2] =~ m!^(\S+):([^\s:]+):([+-]):(\d+):(\d+)$! or die "bad ID [$ta[2]]\n"; 
	my ($mrnaID, $seqID, $str, $pS, $pE) = ($1, $2, $3, $4, $5); 
	my ($seqS, $seqE, $seqStr); 
	if ($str eq '+') {
		$seqS = $pS + $ta[3] - 1; 
		$seqE = $pS + $ta[4] - 1; 
		$seqStr = $str; 
	} elsif ($str eq '-') {
		# $seqS = $pS + ($pE-$pS+1-$ta[3]); 
		# $seqE = $pS + ($pE-$pS+1-$ta[4]); 
		$seqE = $pE - $ta[3] + 1; 
		$seqS = $pE - $ta[4] + 1; 
		$seqStr = $str; 
	} else {
		die "$_\n"; 
	}
	if ($ta[5] eq '-') {
		$seqStr =~ tr/+-/-+/; 
	}
	print join("\t", @ta, $mrnaID, $seqID, $seqS, $seqE, $seqStr)."\n"; 
}
# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
# yytNOR		MELO3C024429T1:Cmel351_Chr01:+:35115704:35118203	349	359	+	16.7447	8.7e-08	1	ACACGTCACCT
# yytNOR		MELO3C000679T1:Cmel351_Chr00:-:17643862:17646361	429	439	-	16.7447	8.7e-08	1	ACACGTCACCT
