#!/usr/bin/perl 
use strict; 

my $trim_tail = 70; 
my $min_scf_len = 500; 
my $max_gap = 500; 
my $min_tail_len = 300; 

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
SEQ: 
for (keys %seq) {
	$_ =~ m/^DropGap/ or next; 
	$seq{$_} =~ s/\s//g; 
	while ( $seq{$_} =~ s/^[ATCGatgc]{0,$trim_tail}[nN]+//o ) { 1; }
	while ( $seq{$_} =~ s/[nN]+[ATGCatgc]{0,$trim_tail}$//o ) { 1; }
	length($seq{$_}) < $min_scf_len and next; 
	my @a1 = split(/[nN]{$max_gap,}/, $seq{$_}); 
	for my $ss (@a1) {
		while ( $ss =~ s/^[ATCGatgc]{0,$trim_tail}[nN]+//o ) { 1; } 
		while ( $ss =~ s/[nN]+[ATGCatgc]{0,$trim_tail}$//o ) { 1; }
		length($ss) < $min_scf_len and next; 
		if ( $ss =~ m/^([ATGCatgc]+)[nN]+([ATGCatgc]+)$/ ) {
			my ($seqH, $seqT) = ($1, $2); 
			if ( length($seqH) < $min_tail_len or length($seqT) < $min_tail_len) {
				for my $s1 ($seqH, $seqT) {
					length( $s1 ) < $min_scf_len and next; 
					$sc_ct ++; 
					my $tseq = $s1; 
					$tseq =~ s/(.{100})/$1\n/g; chomp($tseq); 
					print ">dropGap.$sc_ct [$_]\n$tseq\n"; 
				}
				next SEQ; 
			}
		}
		
		$sc_ct ++; 
		my $tseq = $ss; 
		$tseq =~ s/(.{100})/$1\n/g; chomp($tseq); 
		print ">dropGap.$sc_ct [$_]\n$tseq\n"; 
	}
}

