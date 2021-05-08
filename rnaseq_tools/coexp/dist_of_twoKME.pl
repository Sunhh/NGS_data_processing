#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
!@ARGV and die "perl $0 tag1 file1_KME.txt tag2 file2_KME.txt\n"; 

my $tag1 = shift; 
my $f1   = shift; 
my $tag2 = shift; 
my $f2   = shift; 

my @d1 = &fileSunhh::load_tabFile($f1, 0); 
my @d2 = &fileSunhh::load_tabFile($f2, 0); 

my %simMat; 
my (%g1, %g2); 
for (@d1) {
	$_->[2] eq 'ModuleID' and next; 
	$g1{ $_->[2] }{ $_->[0] } = $_->[1]; 
}
for (@d2) {
	$_->[2] eq 'ModuleID' and next; 
	$g2{ $_->[2] }{ $_->[0] } = $_->[1]; 
}
my @mod1 = sort keys %g1; 
my @mod2 = sort keys %g2; 
for my $m1 (@mod1) {
	my @g1_eleID = sort keys %{$g1{$m1}}; 
	for my $m2 (@mod2) {
		my @g2_eleID = keys %{$g2{$m2}}; 
		my $totalV = 0; 
		my $simV   = 0; 
		for my $e2 (@g2_eleID) {
			defined $g1{$m1}{$e2} and $simV ++; 
			$totalV ++; 
		}
		for my $e1 (@g1_eleID) {
			defined $g2{$m2}{$e1} and next; 
			$totalV ++; 
		}
		$simMat{$m1}{$m2} = $simV/$totalV; 
	}
}
print STDOUT join("\t", 'MM', (map { "${tag2}_$_" } @mod2), (map { "${tag1}_$_" } @mod1))."\n"; 
for my $m2 (@mod2) {
	my @ss; 
	for my $m2_1 (@mod2) {
		my $k=0; 
		$m2 eq $m2_1 and $k = 1; 
		push(@ss, $k); 
	}
	print STDOUT join("\t", "${tag2}_$m2", @ss, (map { $simMat{$_}{$m2} } @mod1))."\n"; 
}
for my $m1 (@mod1) {
	my @ss; 
	for my $m1_1 (@mod1) {
		my $k = 0; 
		$m1 eq $m1_1 and $k = 1; 
		push(@ss, $k); 
	}
	print STDOUT join("\t", "${tag1}_$m1", @{$simMat{$m1}}{@mod2}, @ss)."\n"; 
}


