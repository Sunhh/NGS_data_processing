#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 scrn.run_struct | less -S\n"; 

my @aa; 
while (<>) {
	m!CMD! or next; 
	chomp; 
	s!structure!\tstructure!; 
	s!_(\d+)$!_$1\t$1!; 
	my @ta = split(/\t/, $_); 
	defined $ta[2] or die "bad line: $_\n"; 
	push(@aa, [@ta]); 
}
for (sort { $a->[2] <=> $b->[2] || $a->[2] cmp $b->[2] } @aa){
	print join("\t", @$_)."\n"; 
}

