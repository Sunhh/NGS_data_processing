#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 r6_maker_final.prot.b2g.annot > r6_maker_final.prot.b2g.annot.jnTbl\n"; 

my %gene_infor; 
my %cnt; 
while (<>) {
	m/^\s*(#|$)/ and next; 
	$cnt{'line'}++; 
	chomp; 
	my @ta = split(/\t/, "$_\n"); 
	chomp($ta[-1]); 

	$cnt{'gene_order'}{$ta[0]} //= $cnt{'line'}; 
	if (defined $ta[2]) {
		$gene_infor{$ta[0]}{'def'}{$ta[2]} //= [ $cnt{'line'}, $ta[2] ]; 
	}
	if ($ta[1] =~ m/^GO:\d+$/) {
		$gene_infor{$ta[0]}{'go' }{$ta[1]} //= [ $cnt{'line'}, "$ta[1] ()" ]; 
	} elsif ($ta[1] =~ m/^EC:[\d.]+$/) {
		$gene_infor{$ta[0]}{'ec' }{$ta[1]} //= [ $cnt{'line'}, "$ta[1] ()" ]; 
	} else {
		&stopErr("[Err] Unknown string [$ta[1]] in line: $_\n"); 
	}
}

print STDOUT join("\t", qw/GeneID Definition GOs ECs/)."\n"; 
for my $gid (sort keys %{$cnt{'gene_order'}}) {
	for my $tk (qw/def go ec/) {
		$gene_infor{$gid}{$tk} //= { '-fake' => [-1, ""] }; 
	}
	print STDOUT join("\t", 
	  $gid, 
	  join(';; ', map { $_->[1] } sort { $a->[0] <=> $b->[0] } values %{$gene_infor{$gid}{'def'}}), 
	  join(', ',  map { $_->[1] } sort { $a->[0] <=> $b->[0] } values %{$gene_infor{$gid}{'go' }}), 
	  join(', ',  map { $_->[1] } sort { $a->[0] <=> $b->[0] } values %{$gene_infor{$gid}{'ec' }}), 
	)."\n"; 
}

