#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"shrt_len:i", # 9 
	"shrt_col:i", # 0 
); 

$opts{'shrt_len'} //= 9; 
$opts{'shrt_col'} //= 0; 

my $help_txt = <<HH; 

perl $0 long_indvID_table

-shrt_len           [$opts{'shrt_len'}]
-shrt_col           [$opts{'shrt_col'}]

I will shorten the first column to [$opts{'shrt_len'}] characters. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %h; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[ $opts{'shrt_col'} ] = substr($ta[ $opts{'shrt_col'} ], 0, $opts{'shrt_len'}); 
	my $tk = $ta[ $opts{'shrt_col'} ]; 
	my $suff = "a"; 
	while (defined $h{$tk}) {
		$suff++; 
		$tk = "$ta[ $opts{'shrt_col'} ]$suff"; 
	}
	$h{$tk} = 1; 
	$ta[ $opts{'shrt_col'} ] = $tk; 
	print join("\t", @ta)."\n"; 
}
