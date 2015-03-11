#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"db:s", 
	"inFq1:s@", 
	"inFq2:s@", 
	"outFile:s@", 
	"para_bwt:s", # Default ''
	"exe_bowtie:s", # bowtie
	"exe_bowtieBuild:s", # bowtie-build
	"doBuild!", 
); 

sub usage {
	print <<HH; 
######################################################################
#  perl $0 -db bowtie_database -inFq1 input1.fq -outFile out1 -inFq1 input2.fq -outFile out2
#
#  -inFq2 
#  -para_bwt      [''] Something to consider: '-p cpuN -S -v 0'
#
#  -exe_bowtie       [bowtie]
#  -exe_bowtieBuild  [bowtie-build]
#  -doBuild          [FALSE] Do bowtie-build if given. 
######################################################################
HH
	exit; 
}

$opts{'help'} and &usage(); 

$opts{'exe_bowtie'} //= 'bowtie'; 
$opts{'exe_bowtieBuild'} //= 'bowtie-build'; 
$opts{'para_bwt'} //= ''; 

my $db = $opts{'db'} // &stopErr("[Err] -db = ?\n");
$opts{'doBuild'} and &exeCmd("$opts{'exe_bowtieBuild'} $db $db"); 

defined $opts{'inFq1'} or &stopErr("[Err] -inFq1 = ?\n"); 
my @inFq1 = @{$opts{'inFq1'}}; 
my @inFq2; 
my @outF; 
defined $opts{'inFq2'} and @inFq2 = @{$opts{'inFq2'}}; 
defined $opts{'outFile'} and @outF = @{$opts{'outFile'}}; 
for (my $i=0; $i<@inFq1; $i++) {
	defined $inFq2[$i] or $inFq2[$i] = ''; 
	( defined $inFq1[$i] and $inFq1[$i] ne '' ) or do { &tsmsg("[Err] Skip -inFq1=[$inFq1[$i]]\n"); next; }; 
	defined $outF[$i] or $outF[$i] = ''; 
	if ( $inFq2[$i] eq '' ) {
		&exeCmd_1cmd("$opts{'exe_bowtie'} $opts{'para_bwt'} $db $inFq1[$i] $outF[$i]"); 
	} else {
		&exeCmd_1cmd("$opts{'exe_bowtie'} $opts{'para_bwt'} $db -1 $inFq1[$i] -2 $inFq2[$i] $outF[$i]"); 
	}
}

