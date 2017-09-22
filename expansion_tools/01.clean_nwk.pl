#!/usr/bin/perl
use strict; 
use warnings; 
use Text::Balanced qw( extract_bracketed ); 
use LogInforSunhh; 

# http://www.perlmonks.org/?node_id=547596

while (<>) {
	chomp; 
	my $n = &rmStatNum($_); 
	my $m = &fmtBranch($n); 
	print "$m\n"; 
}

# (((Sly:0.11746590,Vvi:0.07461103)0.7300:0.02307856,Ath:0.14512370)1.0000:0.03408003,(SpiOl:0.05289341,Bvu:0.04303894)1.0000:0.05605601);

sub rmStatNum {
	while ($_[0] =~ s/\)\s*[\d.]+/)/g) {
	# while ($_[0] =~ s/:[\d.]+//g or $_[0] =~ s/\)\s*[\d.]+/)/g) {
	}
	return $_[0]; 
}
sub fmtBranch {
	while ($_[0] =~ s/:(\d+\.\d+)/":" . int($1*100)/eg) {
	}
	return $_[0]; 
}
