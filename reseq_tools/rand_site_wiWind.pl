#!/usr/bin/perl
# 20160805 Output lines with raw order. 
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"snp_tbl:s", # Acc131_mask.snp 
	"wind_len:i", # [10000] window length to select SNP; 10000
	"wind_num:i", # [1] SNP number to be selected from a window. 
	"noHeader!",  
); 

$opts{'wind_len'} //= 10000; 
$opts{'wind_num'} //= 1; 

my $help_txt = <<HH; 

perl $0 -snp_tbl Acc131_mask.snp   

-wind_len       [$opts{'wind_len'}]
-wind_num       [$opts{'wind_num'}]

-noHeader       [Boolean] Provide this parameter if there isn't header in snp table. 

HH

defined $opts{'snp_tbl'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $fh_snpTbl = &openFH( $opts{'snp_tbl'}, '<' ); 
my %stored; 
my %cID_order; 
my $header; 
$opts{'noHeader'} or do { $header = readline($fh_snpTbl); print STDOUT $header;  }; 
my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>1e6 ); 
while (readline($fh_snpTbl)) {
	&fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Reading $. line.\n"); 
	chomp; 
	m!^(\S+)\t(\d+)(?:\t|$)! or die "$_\n"; 
	my ($cID, $cPos) = ($1, $2); 
	my $pIdx = int( $cPos / $opts{'wind_len'} ); 
	push(@{$stored{$cID}[$pIdx]}, [$_, $.]); 
	$cID_order{$cID} //= $.; 
}
close($fh_snpTbl); 

&tsmsg("[Rec] Randomly selecting.\n"); 
for my $cID ( sort { $cID_order{$a} <=> $cID_order{$b} } keys %stored ) {
	&tsmsg("[Msg]   Processing [$cID]\n"); 
	for (my $i=0; $i<@{$stored{$cID}}; $i++) {
		( defined $stored{$cID}[$i] and @{ $stored{$cID}[$i] } > 0 ) or next; 
		my @tb = @{ $stored{$cID}[$i] }; 
		my @slct; 
		my $t_cnt = 0; 
		while ( @tb > 0 ) {
			my $j = int( rand($#tb+1) ); 
			push(@slct, splice( @tb, $j, 1 )); 
			$t_cnt ++; 
			$t_cnt >= $opts{'wind_num'} and last; 
		}
		@slct > 0 or &stopErr("[Err] No sites selected for window [$cID - $i]\n"); 
		for my $tr ( sort { $a->[1] <=> $b->[1] } @slct ) {
			print STDOUT "$tr->[0]\n"; 
		}
	}
}


&tsmsg("[Rec] All done [$0]\n"); 

