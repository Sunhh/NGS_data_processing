#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"digits:i", # 2
); 

-t and !@ARGV and die "\nperl $0 P1g_perGene.sense.cnt.3_noSum_keep_wiSizeFact\n\n"; 

$opts{'digits'} //= 2; 

my @sf; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		print STDOUT "$_\n"; 
		next; 
	}
	if ($. == 2) {
		$ta[0] =~ m/sizeFactor/i or die "The second line should be sizeFactor\n"; 
		@sf = @ta; 
		next; 
	}
	for (my $i=1; $i<@ta; $i++) {
		$ta[$i] = sprintf("%0.$opts{'digits'}f", $ta[$i]/$sf[$i]); 
	}
	print STDOUT join("\t", @ta)."\n"; 
}

