#!/usr/bin/perl
use strict; 
use warnings; 
# Add new ID to .loc table. 
# Format of .loc : ID \\t Start \\t End \\n
# New format :   : ID \\t Start \\t End \\t Pref_ID_Start_End \\n 

!@ARGV and die "perl $0 pref in_raw.loc\n"; 

my $pref = shift; 

my %used; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $new_id = "${pref}_$ta[0]_$ta[1]_$ta[2]"; 
	defined $used{$new_id} and die "repeat new id=$new_id\n"; 
	$used{$new_id} = 1; 
	print STDOUT join("\t", @ta[0,1,2], $new_id)."\n"; 
}

