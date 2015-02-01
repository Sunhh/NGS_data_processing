#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"path_repClass:s", "para_RC:s", 
	"inFa:s", "tagUnknown:s", 
	"noChopID!", 
	"out:s", 
); 

sub usage {
	print STDOUT <<HH; 
################################################################################
# perl $0 -inFa toBeClassified.fas -tagUnknown LTR
#
# -help
# -path_repClass [RepeatClassifier]
# -para_RC       [-engine ncbi]
#
# -out           [inFa.classified]
#
# -noChopID      [boolean] Do not chop ID if given. 
#
################################################################################
HH
	exit 1; 
}
$opts{'path_repClass'} = $opts{'path_repClass'} // '/share/app/Annotation/repeatmodeler/RepeatModeler/RepeatClassifier'; 
$opts{'para_RC'} = $opts{'para_RC'} // '-engine ncbi'; 
$opts{'tagUnknown'} = $opts{'tagUnknown'} // 'LTR'; 

$opts{'help'} and &usage(); 
defined $opts{'inFa'} and -e $opts{'inFa'} or &usage(); 
$opts{'out'} = $opts{'out'} // "$opts{'out'}.classified"; 


my $tmp_file = 'tmp_toClass.fa'; 
unlink($tmp_file); 

unless ( $opts{'noChopID'} ) {
	open IN,'<',"$opts{'inFa'}" or &stopErr("[Err] $!\n"); 
	open OO,'>',"$tmp_file" or &stopErr("[Err] $1\n"); 
	while (<IN>) {
		chomp; 
		if (m/^>/) {
			m/^>(\S+)/ or &stopErr("[Err] Bad format of $_\n"); 
			$_ = $1; 
			s!\#.*$!!; 
		}
		print OO "$_\n"; 
	}
	close OO; 
	close IN; 
} else {
	&exeCmd("cat $opts{'inFa'} > $tmp_file"); 
}

&exeCmd("$opts{'path_repClass'} -consensi $tmp_file $opts{'para_RC'}"); 

open FF,'<',"${tmp_file}.classified" or &stopErr("[Err] $!\n"); 
open BB,'>',"$opts{'out'}" or &stopErr("[Err] $!\n"); 
while (<FF>) {
	chomp; 
	if (m/^>/) {
		s/^>(\S+)// or &stopErr("[Err] Bad format of $_\n"); 
		my $h1 = $1; 
		my $tag1 = ''; 
		if ( $h1 =~ s!\#(\S+)$!! ) {
			$tag1 = $1; 
			$tag1 eq 'Unknown' and $tag1 = $opts{'tagUnknown'}; 
		} else {
			&tsmsg("[Wrn] No tag information for ID [$h1]\n"); 
			$tag1 = $opts{'tagUnknown'}; 
		}
		$_ = ">${h1}#${tag1}$_"; 
	}
	print BB "$_\n"; 
}
close BB; 
close FF; 


