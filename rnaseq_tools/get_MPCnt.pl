#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"digits:i", # 0
	"noSizeFactorLine!", 
); 
$opts{'digits'} //= 0; 

my %h; 
my %cnt; 

{
	my $h = <>; 
	chomp($h); 
	my @ta = split(/\t/, $h); 
	my %done; 
	for (my $i=1; $i<@ta; $i++) {
		# $ta[$i] =~ m/^(P[13]g)_total_P[13]_(\S+)_(rep\d+)$/ or die "$ta[$i]\n"; 
		$ta[$i] =~ m/^(P[13]g)_total_(S(?:FL|FR|LV|RT|SD|ST))(F1|P1|P3)_(rep\d+)$/ or next; 
		my ($pT, $tissue, $indv, $rep) = ($1, $2, $3, $4); 
		my $sid = "$tissue\t$rep"; 
		$h{'have'}{$sid}{$indv} = $i; 
		$h{'idx2sid'}{$i} = $sid; 
		$h{'idx2tissue'}{$i} = $tissue; 
		$cnt{$sid} //= 0; 
		( defined $h{'have'}{$sid}{'P1'} and defined $h{'have'}{$sid}{'P3'} ) or next; 
		defined $done{$sid} and next; 
		$done{$sid} = 1; 
		push(@{$h{'paired_sid'}}, [$h{'have'}{$sid}{'P1'}, $h{'have'}{$sid}{'P3'}, "${pT}_MP_${tissue}_${rep}"]); 
	}
	@{$h{'paired_sid'}} = sort { $a->[2] cmp $b->[2] } @{$h{'paired_sid'}}; 
	print STDOUT join( "\t", "Gene_ID", map { $_->[2] } @{$h{'paired_sid'}} ). "\n"; 
	$opts{'noSizeFactorLine'} or print STDOUT join( "\t", "sizeFactor", map { 1 } @{$h{'paired_sid'}} ). "\n"; 
}

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	my @out = ($ta[0]); 
	for (my $i=0; $i < @{$h{'paired_sid'}}; $i++ ) {
		my $t1 = $h{'paired_sid'}[$i]; 
		my ($i1, $i2) = @$t1; 
		my $md = sprintf("%0.$opts{'digits'}f", ($ta[$i1]+$ta[$i2])/2); 
		push(@out, $md); 
	}
	print join("\t", @out)."\n"; 
}
