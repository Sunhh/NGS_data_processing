#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
use LogInforSunhh; 
my %opts; 
GetOptions(\%opts, 
	"bigCol:i", 
	"smallCol:i", 
	"ratio:f", 
	"help!", 
); 

$opts{'ratio'} //= 1; 

my $help_txt = <<HH; 

perl $0    cafe_input_grp01.report.cafe.sol.tab   -ratio 1   -bigCol 3 -smallCol 4  > cafe_input_grp01.report.cafe.sol.tab.sol_big

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

# [Sunhh@whale run03_Nov_24]$ head cafe_input_grp01.report.cafe.sol.tab
# Description     ID      osa     sol     bvu     vvi
# ORTHOMCL4       ORTHOMCL4       22      49      20      76
# ORTHOMCL8       ORTHOMCL8       31      25      2       40
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[0] eq 'Description' or $ta[2] !~ m/^\d+$/) {
		&tsmsg("[Wrn] Skip line: $_\n"); 
		print STDOUT "$_\n"; 
		next; 
	}
	my $is_good = 0; 
	if ($opts{'ratio'} == 1) {
		$ta[$opts{'bigCol'}] > $ta[$opts{'smallCol'}] and $is_good = 1; 
	} elsif ($opts{'ratio'} > 1) {
		$ta[$opts{'bigCol'}] > 0 and $ta[$opts{'bigCol'}] >= $ta[$opts{'smallCol'}] * $opts{'ratio'} and $is_good = 1; 
	} else {
		&stopErr("[Err] bad value of -ratio $opts{'ratio'}\n"); 
	}
	$is_good == 1 and print STDOUT "$_\n"; 
}

