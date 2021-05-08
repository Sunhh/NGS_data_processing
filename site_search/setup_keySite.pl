#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 fimoOut_allTrans/fimo.tsv.info > fimoOut_allTrans/fimo.tsv.info.keysite\n"; 

print join("\t", qw/chr pos good MotifStart MotifEnd MotifStr PosMotif mrnaID/)."\n"; 
while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	$ta[0] eq 'motif_id' and next; 
	my ($seqID, $pS, $pE, $pStr, $mrnaID) = @ta[ 11, 12, 13, 14, 10 ]; 
	if ($pStr eq '+') {
		print join("\t", $seqID, $pS+3-1, "A",   $pS, $pE, $pStr, 3, $mrnaID)."\n"; 
		print join("\t", $seqID, $pS+4-1, "C",   $pS, $pE, $pStr, 4, $mrnaID)."\n"; 
		print join("\t", $seqID, $pS+5-1, "G",   $pS, $pE, $pStr, 5, $mrnaID)."\n"; 
		print join("\t", $seqID, $pS+8-1, "A,T", $pS, $pE, $pStr, 8, $mrnaID)."\n"; 
	} elsif ($pStr eq '-') {
		print join("\t", $seqID, $pE-3+1, "T",   $pS, $pE, $pStr, 3, $mrnaID)."\n"; 
		print join("\t", $seqID, $pE-4+1, "G",   $pS, $pE, $pStr, 4, $mrnaID)."\n"; 
		print join("\t", $seqID, $pE-5+1, "C",   $pS, $pE, $pStr, 5, $mrnaID)."\n"; 
		print join("\t", $seqID, $pE-8+1, "T,A", $pS, $pE, $pStr, 8, $mrnaID)."\n"; 
	} else {
		die "$_\n"; 
	}
}


# 0 motif_id	yytNOR
# 1 motif_alt_id	
# 2 sequence_name	MELO3C024429T1:Cmel351_Chr01:+:35115704:35118203
# 3 start	349
# 4 stop	359
# 5 strand	+
# 6 score	16.7447
# 7 p-value	8.7e-08
# 8 q-value	1
# 9 matched_sequence	ACACGTCACCT
# 10 mrnaID	MELO3C024429T1
# 11 SeqID	Cmel351_Chr01
# 12 PromoterStart	35116052
# 13 PromoterEnd	35116062
# 14 SeqStr	+
