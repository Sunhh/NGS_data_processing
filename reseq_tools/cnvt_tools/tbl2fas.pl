#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
); 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 

$opts{'startColN'} //= 2; 

my $help_txt = <<HH ; 

perl $0 in_snp.tbl > in_snp.fasta

Need a headere line. 

-help
-startColN       [$opts{'startColN'}]

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $hl = <>; 
chomp($hl); 
my @hh = split(/\t/, $hl); 
my @seq; 

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	for (my $i=$opts{'startColN'}; $i<@ta; $i++) {
		if ( $ta[$i] eq '*' ) {
			$ta[$i] = '-'; 
		} else {
			$ta[$i] = $st_obj->SingleChar($ta[$i], 'maxAlleleN'=>2); 
		}
		$seq[$i] .= $ta[$i]; 
	}
}

for (my $i=$opts{'startColN'}; $i<@hh; $i++) {
	print STDOUT ">$hh[$i]\n"; 
	$seq[$i] =~ s/\s//g; 
	$seq[$i] =~ s/(\S{100})/$1\n/g; 
	chomp($seq[$i]); 
	print STDOUT "$seq[$i]\n"; 
}



