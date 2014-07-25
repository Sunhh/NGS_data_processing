#!/usr/bin/perl 
use strict; 

my $min_len = 150; 
my $max_gap = 700; 

my (%seq, $kk); 
while (<>) {
	if (m/^\s*\>/) {
		m/^\s*\>(\S+)/ or die "$_\n"; 
		$kk = $1; 
	} else {
		$seq{$kk} .= $_; 
	}
}

my $sc_ct = 0; 
for (keys %seq) {
	$seq{$_} =~ s/\s//g; 
	$seq{$_} =~ s/^[nN]+//; 
	$seq{$_} =~ s/[nN]+$//; 
	length($seq{$_}) < $min_len and next; 
	my @a1 = split(/[nN]{$max_gap,}/, $seq{$_}); 
	for my $ss (@a1) {
		length($ss) < $min_len and next; 
		$sc_ct ++; 
		my $tseq = $ss; 
		$tseq =~ s/(.{100})/$1\n/g; chomp($tseq); 
		print ">$sc_ct [$_]\n$tseq\n"; 
	}
}

