#!/usr/bin/perl
use strict; 
use warnings; 

my $avgC1 = 5; # high
my $avgC2 = 6; # low

-t and !@ARGV and die "paste filt1_w50ks5k.pi.gCC.PIavg filt1_w50ks5k.pi.gCA.PIavg | deal_table.pl -column 0-5,11 | perl $0 > CC_CA.avgComp\n"; 

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	if ($ta[0] eq 'ChromID') {
		print STDOUT join("\t", qw/CHROM BIN_START BIN_END BpCnt perKb_High perKb_Low MEAN_Est/)."\n"; 
		next; 
	}
	my $est = ( $ta[$avgC1]+$ta[$avgC2] > 0 ) ? ($ta[$avgC1]/($ta[$avgC1]+$ta[$avgC2])) : 'NA' ; 
	print STDOUT join("\t", @ta[0,1,2,4,$avgC1, $avgC2], $est)."\n"; 
}
