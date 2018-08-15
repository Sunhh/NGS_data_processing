#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 

my $hh = <<"H1"; 

perl $0 S636_fixNH.bam.cntSens S636_fixNH.bam.cntAnti S636_fixNH.bam_pref > S636_fixNH.bam.SensByAnti 

H1

@ARGV >= 3 or &LogInforSunhh::usage($hh); 

my $fn1 = shift; 
my $fn2 = shift; 
my $pref = shift; 

my ($sum1, $sum2); 

for (map { $_->[1] } grep { !($_->[0] =~ m!^__!)} &fileSunhh::load_tabFile($fn1)) {
	$sum1 += $_; 
}
for (map { $_->[1] } grep { !($_->[0] =~ m!^__!)} &fileSunhh::load_tabFile($fn2)) {
	$sum2 += $_; 
}

print STDOUT join("\t", qw/inPref sumSens sumAnti SensByAnti/)."\n"; 
my $r = ($sum2 > 0) ? sprintf("%0.4f", $sum1/$sum2) : -1 ; 
print STDOUT join("\t", $pref, $sum1, $sum2, $r)."\n"; 

