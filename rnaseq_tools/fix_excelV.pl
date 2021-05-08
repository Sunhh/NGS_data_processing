#!/user/bin/perl
use strict; 
use warnings; 
use fileSunhh; 

# Min value that I want to use in Excel is '9E-308'; 
-t and !@ARGV and die "perl $0 compList_inFruit_01.FDR > compList_inFruit_01.FDR_fixSmallV\n"; 

while (<>) {
	chomp; 
	my @ta = &splitL("\t", $_); 
	for my $tb (@ta) {
		if ($tb =~ m!^([+-]?)(\d+(?:\.\d+)?)e\-(\d+)$!i) {
			my ($str, $prevE, $afterE) = ($1, $2, $3); 
			if      ( $afterE > 308 ) {
				$tb = "9E-308";
			} elsif ( $afterE == 308 and $prevE < 9 ) {
				$tb = "9E-308"; 
			} else {
			}
		}
	}
	print STDOUT join("\t", @ta)."\n"; 
}

