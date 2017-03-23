#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 rm_list in_gff\n"; 

my $lisF =shift; 
my $gffF = shift; 
open LF,'<',"$lisF" or die; 
my %rmid; 
while (<LF>) {
	chomp; 
	m/^\s*(#|$)/ and next; 
	my @ta = split(/\t/, $_); 
	$rmid{$ta[0]} = 1; 
}
close LF; 
open GF,'<',"$gffF" or die; 
while (<GF>) {
	chomp; 
	if ( m/^\s*(#|$)/ ) {
		print "$_\n"; 
		next; 
	} 
	my @ta = split(/\t/, $_); 
	if ($ta[8] =~ m/(?:^|;|\s)ID=([^\s;]+)/i) {
		defined $rmid{$1} and next; 
	}
	if ($ta[8] =~ m/(?:^|;|\s)Parent=([^\s;]+)/i) {
		defined $rmid{$1} and next; 
	}
	print "$_\n"; 
}
close GF; 
