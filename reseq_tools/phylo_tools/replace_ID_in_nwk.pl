#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

my $help_txt = <<HH; 

perl $0 old2new_ID in.nwk > new.nwk 

HH

@ARGV == 2 or &LogInforSunhh::usage($help_txt); 

my %id2new = %{ &load_id2new($ARGV[0]) }; 
open AA,'<',"$ARGV[1]" or die; 
while (<AA>) {
	chomp; 
	my @posList; 
	pos($_) = 0; 
	while ($_ =~ m/\G(?:.*?[\(,])([^\s\(\),:]+)(?:\:)/gs) {
		push( @posList, [ $-[1]+1, $+[1], $1 ] ); 
		pos($_) = $+[1]; 
	}
	my $prevE = 0; 
	my @ele; 
	for my $ar (@posList) {
		my ($cS, $cE, $match) = @$ar; 
		if ($cS-1 > $prevE) {
			push( @ele, substr($_, $prevE, $cS-1-$prevE) ); 
		}
		my $toChg = substr($_, $cS-1, $cE-$cS+1); 
		defined $id2new{$toChg} and $toChg = $id2new{$toChg}; 
		push( @ele, $toChg ); 
		$prevE = $cE; 
	}
	if ( length($_) >= $prevE+1 ) {
		push( @ele, substr($_, $prevE, length($_)-$prevE) ); 
	}
	print join('', @ele);
}
close AA; 

sub load_id2new {
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$back{$ta[0]} = $ta[1]; 
	}
	close($fh); 
	return(\%back); 
}
