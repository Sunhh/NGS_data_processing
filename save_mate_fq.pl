#!/usr/bin/perl -w 
use strict; 

@ARGV == 2 or @ARGV == 3 or die "perl $0 mate1.fq mate2.fq [5000000]\n"; 

my $f1 = shift; 
my $f2 = shift; 
my $max_num = shift; 
(defined $max_num and $max_num > 1) or $max_num = 5*(10**6); 
$max_num = int($max_num); 
$max_num > 0 or die "max_num >= 1 !\n"; 

if ($f1 =~ /\.gz$/) {
	open F1,"-|", "gzip -cd $f1" or die; 
}else{
	open F1,'<',"$f1" or die; 
}

open O1,'>',"$f1.mate" or die; 
open O2,'>',"$f2.mate" or die; 
my (%record1, $num1, $k1); 
$num1 = 0; 
my $pair_num1 = 0; 
my $pair_num2 = 0; 
while (my $l1 = <F1>) {
	my $n1 = $. % 4; # one sequence in each line for fastq format. 
	if ($n1 == 1) {
		$l1 = "\@$l1"; 
		if ($num1 == $max_num) {
warn "[Stat]Checking pairs at line $. in file1: $f1\n"; 
			if ($f2 =~ /\.gz$/) {
				open F2,"-|","gzip -cd $f2" or die; 
			}else{
				open F2,'<',"$f2" or die; 
			}
			my $is_out2 = 0; 
			while (my $l2 = <F2>) {
$. % 1000000 == 1 and warn "[Stat]Checking line $. in file2: $f2\n"; 
				my $n2 = $. % 4; 
				if ($n2 == 1) {
					$l2 = "\@$l2"; 
					$is_out2 = 0; 
					(my $k2) = ($l2 =~ /^\@(\S+)\/2/) or die "f2:$l2\n"; 
					if (defined $record1{$k2}) {
						$is_out2 = 1; 
						print O1 $record1{$k2}; 
						$pair_num1 ++; 
					}
				}
				$is_out2 == 1 and do { print O2 $l2; $pair_num2 ++; }; 
			}
			close F2; 
			$num1 = 0; 
			undef(%record1); 
		}
		
		$num1 ++; 
		($k1) = ($l1 =~ /^\@(\S+)\/1/) or die "f1:$l1\n"; 
		$record1{$k1} = $l1; 
	}else{
		$record1{$k1} .= $l1; 
	}
}
close F1; 
if ($num1 > 0) {
warn "[Stat]Checking rest lines.\n"; 
	if ($f2 =~ /\.gz$/) {
		open F2,"-|", "gzip -cd $f2" or die; 
	}else{
		open F2,'<',"$f2" or die; 
	}
	my $is_out2 = 0; 
	while (my $l2 = <F2>) {
$. % 1000000 == 1 and warn "[Stat]Checking line $. in file2: $f2\n"; 
		my $n2 = $. % 4; 
		if ($n2 == 1) {
			$l2 = "\@$l2"; 
			$is_out2 = 0; 
			(my $k2) = ($l2 =~ /^\@(\S+)\/2/) or die "f2:$l2\n"; 
			if (defined $record1{$k2}) {
				$is_out2 = 1; 
				print O1 $record1{$k2}; 
				$pair_num1 ++; 
			}
		}
		$is_out2 == 1 and do { print O2 $l2; $pair_num2++; }; 
	}
	close F2; 
	undef(%record1)
}
close O1; 
close O2; 

warn "Game passed!\n"; 


