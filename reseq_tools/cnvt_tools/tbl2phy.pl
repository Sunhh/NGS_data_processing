#!/usr/bin/perl
use strict; 
use warnings; 

my @header; 
my @seqs; 
my $colS = 2; 
while (<>) {
	s/[^\S\t]+$//; 
	my @ta = split(/\t/, $_); 
	if ( $. == 1) {
		@header=@ta; 
		next; 
	}
	for (my $i=$colS; $i<@ta; $i++) {
		$ta[$i] =~ m/^[ATGC]$/ or $ta[$i] = 'N'; 
		length($ta[$i]) == 1 or die "ta[$i]=[$ta[$i]] in line:\n$_\n"; 
		$seqs[$i] .= $ta[$i]; 
	}
}

my $indv_N = $#header-$colS+1; 
my $base_N = length($seqs[$colS]); 

print STDOUT <<OOO; 
   $indv_N   $base_N
OOO

my %used_id; 
for (my $i=$colS; $i<@header; $i++) {
	my $indv_ID = $header[$i]; 
	$indv_ID =~ s!\s!_!g; 
	$indv_ID =~ tr!)(][:;,!_______!; 
	length($indv_ID) > 9 and $indv_ID = substr($indv_ID, 0, 9); 
	$indv_ID = sprintf("%-10s", $indv_ID); 
        defined $used_id{$indv_ID} and die "Repeat ID [$indv_ID]\n"; 
	$used_id{$indv_ID} = 1; 
	$seqs[$i] =~ s!\s!!g; 
	print STDOUT "#${indv_ID}$seqs[$i]\n"; 
}

