#!/usr/bin/perl -w 
use strict; 
-t and !@ARGV and die "perl $0 in.scaffold > out.scaffold_ctg\n"; 

my ($k, $seq) = ('', ''); 
while (<>) {
	if (/^>(\S+)/) {
		my $tk = $1; 
		if ($seq !~ /^\s*$/) {
			$seq =~ s/\s//g; 
			my $num = 0; 
			# Method 1. 
			#$seq =~ s/^[Nn]+//; 
			#while ($seq =~ s/^([^Nn]+)[Nn]+//) {
			#	$num ++; 
			#	print STDOUT ">${k}_$num\n$1\n"; 
			#}
			#if ($seq ne '') {
			#	$num ++; 
			#	print STDOUT ">${k}_$num\n$seq\n"; 
			#	$seq = ''; 
			#}
			# Method 2. Should be better. 
			pos($seq) = 0; 
			while ($seq =~ m/\G(?:.*?)([^Nn]+)/gs) {
				$num ++; 
				my ($s, $e, $len) = ($-[1]+1, $+[1], $+[1]-$-[1]); 
				print STDOUT ">${k}_$num [$s,$e,$len]\n$1\n"; 
				pos($seq) = $+[1]; 
			}
		}
		$k = $tk; $seq = ''; 
	}else{
		$seq .= $_; 
	}
}

if ($seq !~ /^\s*$/) {
	$seq =~ s/\s//g; 
	my $num = 0; 
	pos($seq) = 0; 
	while ($seq =~ m/\G(?:.*?)([^Nn]+)/gs) {
		$num ++; 
		my ($s, $e, $len) = ($-[1]+1, $+[1], $+[1]-$-[1]); 
		print STDOUT ">${k}_$num [$s,$e,$len]\n$1\n";
		pos($seq) = $+[1]; 
	}
	$seq = ''; $k = ''; 
}


