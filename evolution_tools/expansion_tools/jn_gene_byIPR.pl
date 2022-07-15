#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "perl $0 in.iprV5.tsv\n"; 

my $fn = shift; 
open F,'<',"$fn" or die; 
my %h; 
while (&wantLineC(\*F)) {
	my @ta = &splitL("\t", $_); 
	( defined $ta[11] and $ta[11] ne '' ) or next; 
	my $tk = "$ta[11] ($ta[12])"; 
	$h{$tk}{$ta[0]} //= $.; 
}
close F; 
print STDOUT join("\t", 'IPRID', $fn)."\n"; 
for my $iprK ( sort keys %h ) {
	my @tb = sort { $a cmp $b } keys %{$h{$iprK}}; 
	print STDOUT join("\t", $iprK, scalar(@tb) . ":" . join(',', @tb))."\n"; 
}

