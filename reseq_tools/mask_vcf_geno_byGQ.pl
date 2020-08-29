#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, "help!", 
	"minGQ:i", 
	"minDP:i", 
); 

$opts{'minGQ'} //= -1; 
$opts{'minDP'} //= -1; 

my $help_txt = <<HH; 
################################################################################
# perl $0 -minGQ 20 lmyPM2019bsa_filtV_PASS.vcf > lmyPM2019bsa_filtV_PASS_mskGQ.vcf
#
# -minGQ      [$opts{'minGQ'}]
#
################################################################################
HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 


my %inf; 
LINE: 
while (<>) {
	chomp; 
	if (m!^\s*#!) {
		m!^\s*##! and do { print STDOUT "$_\n"; next LINE; }; 
		# m!^\s*##! and do { next LINE; }; 
		if ( m!^#CHROM\t! ) {
			print "$_\n"; 
			next LINE; 
		} else {
			&stopErr("[Err] Bad line: $_\n"); 
		}
	}
	my @ta=split(/\t/, $_); 
	my $is_good = 1; 
	my @a_f = split(/:/, $ta[8]); 
	my $j_GT; 
	SAMPLEACC: 
	for (my $i=9; $i<@ta; $i++) {
		my @a_g = split(/:/, $ta[$i]); 
		for (my $j=0; $j<@a_g; $j++) {
			defined $a_f[$j] or last; 
			if ($a_f[$j] eq 'GT') {
				$j_GT //= $j; 
				$a_g[$j] eq './.' and next SAMPLEACC; 
			} elsif ( $a_f[$j] eq 'GQ' ) {
				$a_g[$j] =~ m!^\d+$! or do { &tsmsg("[Wrn] Skip bad line: j=$j:$a_g[$j]: $_\n"); print STDOUT "$_\n";  next LINE; }; 
				if ( $a_g[$j] < $opts{'minGQ'} ) {
					if (!defined $j_GT) {
						for (my $x=0; $x<@a_f; $x++) {
							$a_f[$x] eq 'GT' and do { $j_GT = $x; last; }; 
						}
						defined $j_GT or do { &tsmsg("[Wrn] Skip bad line:noJ: $_\n"); print STDOUT "$_\n";  next LINE; }; 
					}
					$a_g[$j_GT] = './.'; 
				}
				last; 
			}
		}
		$ta[$i] = join(":", @a_g); 
	}
	print STDOUT join("\t", @ta)."\n"; 

}




