#!/usr/bin/perl 
# 2016-12-01 Add -cpuN 
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"cpuN:i", 
	"colN_ref:i", # 2
	"colN_start:i", # 3
	"noSrtAllele!", # Do not sort alleles within genotype by character order. 
); 

$opts{'cpuN'} //= 1; 
$opts{'cpuN'} = int($opts{'cpuN'}); 
$opts{'cpuN'} < 1 and $opts{'cpuN'} = 1; 
$opts{'colN_ref'}   //= 2; 
$opts{'colN_start'} //= 3; 


my $help_txt = <<HH; 

perl $0 in_sMao.snp.cols > in_sMao.snp.tab

-help
-cpuN           [$opts{'cpuN'}] 

-colN_ref       [$opts{'colN_ref'}]
-colN_start     [$opts{'colN_start'}]

-noSrtAllele    [Boolean] If given, I won't sort alleles within genotype, 
                          so the 'AT' will be 'A/T', and 'TA' will be 'T/A'. 
                          If not given, 'AT' or 'TA' will always be 'A/T'. 

Format of in_sMao.snp.cols : 
  chr \\t pos \\t base(ref) \\t sample1 \\t sample2 ... 
  c1  \\t 100 \\t A         \\t AT      \\t G
  c2  \\t 10  \\t G         \\t N       \\t *
  c3  \\t 5   \\t T         \\t A+ATG   \\t *A+ATG
  c4  \\t 8   \\t N         \\t A*      \\t R
  c4  \\t 10  \\t A         \\t AGT     \\t H

Format of vcf.tab : 
  chr \\t pos \\t base(ref) \\t sample1    \\t sample2 ...
  c1  \\t 100 \\t A         \\t A/T        \\t G/G
  c2  \\t 10  \\t G         \\t ./.        \\t */*
  c3  \\t 5   \\t T         \\t AATG/AATG  \\t AATG/*
  c4  \\t 8   \\t N         \\t A/*        \\t A/G
  c4  \\t 10  \\t A         \\t ./.        \\t ./.

By default, I'll keep column 'base' as single letter, and 

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
-t and !@ARGV and &LogInforSunhh::usage($help_txt); 

my @InFp = (); 
if ( !@ARGV ) {
	-t or @InFp = (\*STDIN); 
} else {
	for (@ARGV) {
		push( @InFp, &openFH($_, '<') ); 
	}
}
my $pm; 
$opts{'cpuN'} > 1 and $pm = &LogInforSunhh::get_pm( $opts{'cpuN'} ); 


while (<>) { 
	chomp; 
	my @ta=split(/\t/, $_); 
	if ( $. == 1 ) {
		print join("\t", @ta)."\n"; 
		next; 
	}
	for my $tb (@ta[3..$#ta]) { 
		if      ( $tb =~ s!^([ATGC*])$!$1/$1! ) {
			; 
		} elsif ( $tb =~ s!^[nN]$!./.!i ) {
			; 
		} elsif ( $tb =~ s!^([ATGC*])([ATGC*])$!$1/$2! ) {
			; 
		} elsif ( $tb =~ s!^([ATGC])\+([ATGCN]+)$!$1$2/$1$2!) {
			;
		} elsif ( $tb =~ s!^([ATGC])([ATGC])\+([ATGC]+)$!$1/$2$3! ) {
			; 
		} elsif ( $tb =~ m!^\+!) { 
			$tb = "./."; 
		} else { 
			die "$tb\n"; 
		}
		$tb=join("/", sort ($tb =~ m!^([^/]+)/([^/]+)$!) );  
	} 
	print join("\t", @ta)."\n";  
}

