#!/usr/bin/perl -w 
use strict; 

!@ARGV and die "perl $0 input.fastq.mate input.fastq\n"; 

my $f1 = shift; 
my $f2 = shift; 
if ($f1 =~ /\.gz$/) {
	open F1, "-|", "gzip -cd $f1" or die; 
}else{
	open F1, '<', "$f1" or die; 
}
if ($f2 =~ /\.gz$/) {
	open F2, "-|", "gzip -cd $f2" or die; 
}else{
	open F2, '<', "$f2" or die; 
}
open OS2, '>', "$f2.single" or die; 

my ($k1, $k2); 
$k1 = -1; 
$k2 = -2; 
# my ($record1, $record2); 
my $is_out = 1; 
FILE1: 
while (my $l1 = <F1>) {
	my $n1 = $. % 4; 
$. % 1000000 == 1 and warn "[Stat]$. line in $f1\n"; 
	if ($n1 == 1) {
		# $l1 = "\@$l1"; 
		($k1) = ($l1 =~ /^@(\S+)\/1/); 
		FILE2: 
		while (my $l2 = <F2>) {
			my $n2 = $. % 4; 
			if ($n2 == 1) {
				$l2 = "\@$l2"; 
				($k2) = ($l2 =~ /^@(\S+)\/[12]/); 
				if ($k1 eq $k2) {
					$is_out = 0; 
					last FILE2; 
				}else{
					$is_out = 1; 
				}
			}
			$is_out == 1 and print OS2 $l2; 
		}#file2
	}
}

while (my $l2 = <F2>) {
	my $n2 = $. % 4; 
	if ($n2 == 1) {
		$l2 = "\@$l2"; 
		$is_out = 1; 
	}
	$is_out == 1 and print OS2 $l2; 
}
close F1; 
close F2; 
close OS2; 

