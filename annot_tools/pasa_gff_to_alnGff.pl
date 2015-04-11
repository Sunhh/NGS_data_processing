#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 P1denovoAndGG_pasa.pasa_assemblies.named.gff3 > P1denovoAndGG_pasa.pasa_assemblies.named.fmt.gff3\n"; 

my $inF=shift; 
my %se; 
open F,'<',"$inF" or die; 
while (<F>) {
	chomp; 
	m!^\s*(#|$)! and next; 
	my @ta = split(/\t/, $_); 
	$ta[0] =~ s!\s.*$!!; 
	$ta[8] =~ s!^ID=([^\s;]+);Target=(\S+)!Parent=3:$1;Target=$2! or die "$_\n"; 
	my $id="3:$1"; 
	my $tgt = $2; 
	$se{$id}{tgt} //= $tgt; 
	$se{$id}{tgt} eq $tgt or die "$_\n"; 
	$se{$id}{s} //= $ta[3]; 
	$se{$id}{e} //= $ta[4]; 
	$se{$id}{str} //= $ta[6]; 
	$se{$id}{ln} //= $.; 
	$se{$id}{s} > $ta[3] and $se{$id}{s} = $ta[3]; 
	$se{$id}{e} < $ta[4] and $se{$id}{e} = $ta[4]; 
	$ta[1] = "P1denovoAndGG_pasa"; 
	$ta[2] = "match_part"; 
	$ta[8] =~ s!;$!!; $ta[8] = "$ta[8];"; 
	push(@{$se{$id}{lines}}, [@ta]); 
}
close F; 
for my $id (sort {$se{$a}{ln} <=> $se{$b}{ln}} keys %se) {
	my %th = %{$se{$id}}; 
	print STDOUT join("\t", $th{lines}[0][0], $th{lines}[0][1], "match", $th{s}, $th{e}, '.', $th{str}, '.', "ID=$id;Name=$th{tgt};")."\n"; 
	if ( $th{str} eq '-' ) {
		@{$th{lines}} = sort { $b->[3] <=> $a->[3] } @{$th{lines}}; 
	} elsif ( $th{str} eq '+' ) {
		@{$th{lines}} = sort { $a->[3] <=> $b->[3] } @{$th{lines}}; 
	} else {
		die "str=$th{str}\n"; 
	}
	for my $tr ( @{$th{lines}} ) {
		print STDOUT join("\t", @$tr)."\n"; 
	}
}
