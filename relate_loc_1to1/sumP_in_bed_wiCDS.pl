#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"addCDS!", 
	"locCols:s", # 2,3
); 

!@ARGV and die "perl $0 P1toP3_1to1.maf.rel_loc P1Genom_V1p2.prot_chr.gff3.bed_fit_YZ\n-addCDS       [Boolean]\n-locCols      [2,3]\n\n"; 

my $fn_loc = shift; 
my $fn_bed = shift; 

my ($c_id, $c_p) = (2,3); 
if ( defined $opts{'locCols'} ) {
	$opts{'locCols'} =~ m/^(\d+),(\d+)$/ or die "bad -locCols $opts{'locCols'}"; 
	($c_id, $c_p) = ($1, $2); 
}

&tsmsg("[Msg] Reading $fn_loc\n"); 
open L,'<',"$fn_loc" or die; 
my %cnt; 
$cnt{'cntN_step'} = 5e6; 
my %have; 
while (<L>) {
	&fileSunhh::log_section( $. , \%cnt ) and &tsmsg("[Msg]   Reading $. line.\n"); 
	chomp; m/^\s*$/ and next; 
	my @ta = split(/\t/, $_); 
	$have{$ta[$c_id]}{$ta[$c_p]} = 1; 
}
close L; 
&tsmsg("[Msg] Done $fn_loc\n"); 


my %sum; 
my @arr; 
%cnt=(); 
$cnt{'cntN_step'} = 1000; 
open B,'<',"$fn_bed" or die; 
while (<B>) {
	&fileSunhh::log_section( $. , \%cnt ) and &tsmsg("[Msg]   Reading $. line.\n"); 
	chomp; m/^\s*$/ and next; 
	my @ta = split(/\t/, $_); 
	push(@arr, [@ta]); 
	$sum{$ta[3]} = 0; 
	defined $have{$ta[0]} or next; 
	if ( $opts{'addCDS'} ) {
		for my $tP ( map { m/^(\d+),(\d+)$/ or die "$_\n"; ($1 .. $2) } split(/;/, $ta[6]) ) {
			defined $have{$ta[0]}{$tP} and $sum{$ta[3]} ++; 
		}
	} else {
		for my $tP ( ($ta[1]+1) .. $ta[2] ) {
			defined $have{$ta[0]}{$tP} and $sum{$ta[3]} ++; 
		}
	}
}
close B; 
for (@arr) {
	print join("\t", @$_, $sum{$_->[3]})."\n"; 
}


