#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 in_R1.fastq\n"; 

my $r1pat_set3 = 'ATCTCGTATGCCGTCTTCTGCTTG'; 
my $r1pat_set1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'; 

my $r2pat = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'; 

for my $f (@ARGV) {
	my $cnt_found = 0; 
	my %cnt_n40; 
	my %cnt_n41; 
	my $fh = &openFH($f, '<'); 
	while (my $id = <$fh>) {
		$cnt_found >= 10000 and last; 
		my $seq = <$fh>; 
		<$fh>; <$fh>; 
		if ( $seq =~ m/$r1pat_set3/o ) {
			my $is = 0; 
			if ( $seq =~ m/(.{40}$r1pat_set3)/o ) {
				$cnt_n40{$1} ++; 
				$is = 1; 
			}
			if ( $seq =~ m/(.{41}$r1pat_set3)/o ) {
				$cnt_n41{$1} ++; 
				$is = 1; 
			}
			$is == 1 and $cnt_found ++; 
		}
	}
	my @k_n40 = sort { $cnt_n40{$b} <=> $cnt_n40{$a} } keys %cnt_n40; 
	my @k_n41 = sort { $cnt_n41{$b} <=> $cnt_n41{$a} } keys %cnt_n41; 
	
	print STDERR "prevNbp=40\n"; 
	for (my $i=0; $i<5; $i++) {
		defined $k_n40[$i] or next; 
		print STDERR join("\t", $k_n40[$i], $cnt_n40{ $k_n40[$i] })."\n"; 
	}
	print STDERR "prevNbp=41\n"; 
	for (my $i=0; $i<5; $i++) {
		defined $k_n41[$i] or next; 
		print STDERR join("\t", $k_n41[$i], $cnt_n41{ $k_n41[$i] })."\n"; 
	}
	my $has_n40 = 0; 
	if ( defined $k_n40[1] and $cnt_n40{ $k_n40[0] } > 3 * $cnt_n40{ $k_n40[1] } ) {
		print STDERR "Found R1_pattern with N40 : \t$k_n40[0]\n"; 
		$has_n40 = 1; 
	}
	my $has_n41 = 0; 
	if ( defined $k_n41[1] and $cnt_n41{ $k_n41[0] } > 3 * $cnt_n41{ $k_n41[1] } ) {
		print STDERR "Found R1_pattern with N41 : \t$k_n41[0]\n"; 
		$has_n41 = 1; 
	}
	if ( $has_n40 == 1 and $has_n41 == 0 ) {
		print STDERR "Good R1_pattern : \t$k_n40[0]\n"; 
		print STDOUT join("\t", $f, $k_n40[0])."\n"; 
		if ( $k_n40[0] =~ m!^${r1pat_set1}(\S+)${r1pat_set3}$! ) {
			print STDERR "Good barcode : \t$1\n"; 
		}
	}
}# for my $f (@ARGV) 

