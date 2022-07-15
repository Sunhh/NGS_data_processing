#!/usr/bin/perl
use strict; 
use warnings; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"ref2id:s", 
); 
; 

-t and !@ARGV and die "perl $0 [-ref2id REF] in.vcfTab > out.vcfTab\n"; 


while (<>) {
	chomp; 
	s!^(\S+\t\S+)\t(\S+)\t!! or die "bad line: $_\n"; 
	if ($. == 1) {
		$opts{'ref2id'} //= $2; 
		print STDOUT "$1\t$2\t$opts{'ref2id'}\t$_\n"; 
	} else {
		print STDOUT "$1\t$2\t$2/$2\t$_\n"; 
	}
}

