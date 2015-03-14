#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 in_R1.fastq\n"; 

my $r1pat_set3 = 'ATCTCGTATGCCGTCTTCTGCTTG'; 
my $r1pat_set1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'; 

my $r2pat = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'; 
my $n1_len=40;
my $n2_len=41; 

for my $f (@ARGV) {
	my $cnt_found = 0; 
	my %cnt_n1; 
	my %cnt_n2; 
	my $fh = &openFH($f, '<'); 
	while (my $id = <$fh>) {
		$cnt_found >= 10000 and last; 
		my $seq = <$fh>; 
		<$fh>; <$fh>; 
		if ( $seq =~ m/$r1pat_set3/o ) {
			my $is = 0; 
			if ( $seq =~ m/(.{$n1_len}$r1pat_set3)/o ) {
				$cnt_n1{$1} ++; 
				$is = 1; 
			}
			if ( $seq =~ m/(.{$n2_len}$r1pat_set3)/o ) {
				$cnt_n2{$1} ++; 
				$is = 1; 
			}
			$is == 1 and $cnt_found ++; 
		}
	}
	my @k_n1 = sort { $cnt_n1{$b} <=> $cnt_n1{$a} } keys %cnt_n1; 
	my @k_n2 = sort { $cnt_n2{$b} <=> $cnt_n2{$a} } keys %cnt_n2; 
	
	print STDERR "prevNbp=$n1_len\n"; 
	for (my $i=0; $i<5; $i++) {
		defined $k_n1[$i] or next; 
		print STDERR join("\t", $k_n1[$i], $cnt_n1{ $k_n1[$i] })."\n"; 
	}
	print STDERR "prevNbp=$n2_len\n"; 
	for (my $i=0; $i<5; $i++) {
		defined $k_n2[$i] or next; 
		print STDERR join("\t", $k_n2[$i], $cnt_n2{ $k_n2[$i] })."\n"; 
	}
	my $has_n1 = 0; 
	if ( defined $k_n1[1] and $cnt_n1{ $k_n1[0] } > 3 * $cnt_n1{ $k_n1[1] } ) {
		print STDERR "Found R1_pattern with N$n1_len : \t$k_n1[0]\n"; 
		$has_n1 = 1; 
	}
	my $has_n2 = 0; 
	if ( defined $k_n2[1] and $cnt_n2{ $k_n2[0] } > 3 * $cnt_n2{ $k_n2[1] } ) {
		print STDERR "Found R1_pattern with N$n2_len : \t$k_n2[0]\n"; 
		$has_n2 = 1; 
	}
	if ( $has_n1 == 1 and $has_n2 == 0 ) {
		print STDERR "Good R1_pattern : \t$k_n1[0]\n"; 
		print STDOUT join("\t", $f, $k_n1[0])."\n"; 
		if ( $k_n1[0] =~ m!^${r1pat_set1}(\S+)${r1pat_set3}$! ) {
			print STDERR "Good barcode : \t$1\n"; 
		}
	}
}# for my $f (@ARGV) 

