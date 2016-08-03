#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"showTime:i", 
); 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 

$opts{'startColN'} //= 2; 

my $help_txt = <<HH ; 

perl $0 in_snp.tbl > in_snp.fasta

Need a headere line. 

-help
-startColN       [$opts{'startColN'}]

-showTime        [0] Will report time when 'showTime' number of lines have been processed. 

HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 
$opts{'showTime'} //= 0; 

my $hl = <>; 
chomp($hl); 
my @hh = split(/\t/, $hl); 
my @seq; 

my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>1e5 ); 
$opts{'showTime'} > 0 and $tmp_cnt{'cntN_step'} = $opts{'showTime'}; 
while (<>) {
	$opts{'showTime'} > 0 and &fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading $. line.\n");
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



