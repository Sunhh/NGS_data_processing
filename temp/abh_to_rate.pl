#!/usr/bin/perl
use strict; 
use warnings; 
# [Sunhh@bioinfor01 recomb_rate]$ head geno_list_abh
# C       PM      NEW-2   NEW-4   NEW-5   NEW-7
# 10-1-1  A       A       A       A       A
# 10-1-10 A       H       H       H       H
# 10-1-11 A       H       H       H       H


my %geno; 
my @hh; 
{
	my $l=<>; 
	chomp($l); 
	my @ta=split(/\t/, $l); 
	@hh = @ta; 
}
while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	for (my $i=1; $i<@ta; $i++) {
		$ta[$i] =~ m!^(A|B|H|Da|Rb|\-)$! or die "$ta[$i]\n"; 
		push(@{$geno{$hh[$i]}}, $ta[$i]); 
	}
}
my %calc; 
my @mm = sort keys %geno; 
for (my $i1=0; $i1<@mm; $i1++) {
	for (my $i2=$i1+1; $i2<@mm; $i2++) {
		my ($cnt_ttl, $cnt_cross) = (0,0); 
		for (my $j=0; $j<@{$geno{$mm[$i1]}}; $j++) {
			my $a2 = "$geno{$mm[$i1]}[$j]$geno{$mm[$i2]}[$j]"; 
			if ( $a2 =~ m!\-! ) {
				# Missing
				; 
			} elsif ( $a2 =~ m!^(AA|BB|HH)$! ) {
				# No cross found. 
				$cnt_ttl += 2; 
			} elsif ( $a2 =~ m!^(AH|BH|HA|HB)$! ) {
				$cnt_ttl += 2; 
				$cnt_cross += 1; 
			} elsif ( $a2 =~ m!^(AB|BA)$! ) {
				$cnt_ttl += 2; 
				$cnt_cross += 2; 
			} elsif ( $a2 =~ m!^(DaA|ADa)$! ) {
				# "Da" means dominant phenotype as A; 
				$cnt_ttl += 2; 
				$cnt_cross += (1*2/3); 
			} elsif ( $a2 =~ m!^(DaB|BDa)$! ) {
				$cnt_ttl += 2; 
				$cnt_cross += (2*1/3+1*2/3); 
			} elsif ( $a2 =~ m!^(DaH|HDa)$! ) {
				$cnt_ttl += 2; 
				$cnt_cross += (1*1/3); 
			} elsif ( $a2 =~ m!^(RbA|ARb)$! ) {
				# "Rb" means recessive phenotype as B; 
				$cnt_ttl += 2; 
				$cnt_cross += 2; 
			} elsif ( $a2 =~ m!^(RbB|BRb)$! ) {
				$cnt_ttl += 2; 
			} elsif ( $a2 =~ m!^(RbH|HRb)$! ) {
				$cnt_ttl += 2; 
				$cnt_cross += 1; 
			} else {
				die "unknown [$a2] at [$mm[$i1]-$j][$mm[$i2]-$j]\n"; 
			}
		}
		print STDOUT join("\t", $mm[$i1], $mm[$i2], $cnt_cross/$cnt_ttl*100, $cnt_ttl, $cnt_cross)."\n"; 
	}
}

