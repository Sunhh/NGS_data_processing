#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"perWind!", 
); 
my $help_txt = <<HH; 

perl $0 list_of_fst_perSiteChrPos > merged.ChrPos

-perWind        [Bool]

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %used; 
while (my $fn = <>) {
	chomp($fn); 
	my @ta = split(/\t/, $fn); 
	my $fh = &openFH($ta[0], '<'); 
	while (<$fh>) {
		chomp; 
		my @tb = split(/\t/, $_); 
		my $tk = ($opts{'perWind'}) ? $tb[0] : "$tb[0]\t$tb[1]"; 
		defined $used{$tk} and next; 
		$used{$tk} = 1; 
		print STDOUT "$_\n"; 
	}
	close($fh); 
}

