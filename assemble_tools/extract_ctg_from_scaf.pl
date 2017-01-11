#!/usr/bin/perl -w 
use strict; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"out_agp:s", 
); 
-t and !@ARGV and die "perl $0 in.scaffold > out.scaffold_ctg\n\n-out_agp   [out.scaffold_ctg2SCF.agp]\n"; 

my $ofn_agp = '/dev/null'; 
defined $opts{'out_agp'} and $ofn_agp = $opts{'out_agp'}; 
open OA,'>',"$ofn_agp" or die; 

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
			my $part_number = 0; 
			my (@prev); 
			while ($seq =~ m/\G(?:.*?)([^Nn]+)/gs) {
				$num ++; 
				my ($s, $e, $len) = ($-[1]+1, $+[1], $+[1]-$-[1]); 
				my $ctgID = "${k}_$num"; 
				print STDOUT ">$ctgID [$s,$e,$len]\n$1\n"; 
				$part_number ++; 
				if ( $part_number > 1 ) {
					print OA join("\t", $k, $prev[2]+1, $s-1, $part_number, 'N', 'scaffold', 'yes', 'paired-ends')."\n"; 
					$part_number ++; 
				}
				print OA join("\t", $k, $s, $e, $part_number, 'W', $ctgID, 1, $len, '+')."\n"; 
				@prev = ( $k, $s, $e, $part_number, 'W', $ctgID, 1, $len, '+' ); 
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


