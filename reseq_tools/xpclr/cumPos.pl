#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"out_cumLen:s", 
	"delt_len:i", # 0 
	"joinedID:s", 
	"col_id:i", # 0
	"col_pos:i@", # 1,2
); 

$opts{'delt_len'} //= 0; 
$opts{'joinedID'} //= 'joinedChr'; 
$opts{'col_id'}   //= 0; 
$opts{'col_pos'}  //= [1,2]; 

my $help_txt = <<HH; 
################################################################################
# perl $0 chrLen xpclr_w10ks10k > xpclr_w10ks10k.jnChr
#
# -help
#
# -delt_len       [0] Distance between two adjacent chromosomes. 
# -out_cumLen     [filename] Output a chrLen_cum file for plotting. 
# -joinedID       [$opts{'joinedID'}] 
# 
# -col_id         [$opts{'col_id'}]
# -col_pos        [-col_pos 1 -col_pos 2]
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
!@ARGV and &LogInforSunhh::usage($help_txt); 

my %glob; 
defined $opts{'out_cumLen'} and $glob{'fh_oCumLen'} = &openFH( $opts{'out_cumLen'}, '>' ); 

my $fn_chrCL = shift; 
my $fn_wind  = shift; 
my $new_id = $opts{'joinedID'};  

my %chrCumS = %{ &load_chrCL( $fn_chrCL ) }; 

my $fh = &openFH( $fn_wind, '<' ); 
while (&wantLineC($fh)) {
	my @ta = &splitL("\t", $_); 
	if ( $ta[0] =~ m/^chrID$/i ) {
		print STDOUT join("\t", @ta)."\n"; 
		next; 
	}
	my $rawID = $ta[$opts{'col_id'}]; 
	$rawID =~ m!^(\d+)$! and $rawID = "chr$rawID"; 
	if ( ! defined $chrCumS{$rawID} ) {
		&tsmsg("[Wrn] Skip line : $_\n"); 
		print STDOUT join("\t", @ta)."\n"; 
		next; 
	}
	for my $tb (@ta[@{$opts{'col_pos'}}]) {
		$tb =~ m!^NA$!i and next; 
		$tb = $tb + $chrCumS{$rawID}[0]; 
	}
	$ta[$opts{'col_id'}] = $new_id; 
	print STDOUT join("\t", @ta)."\n"; 
}
close($fh); 


sub load_chrCL {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	my $cum = 0; 
	defined $opts{'out_cumLen'} and print {$glob{'fh_oCumLen'}} join("\t", qw/chrID chrLen chrCumS chrCumE/)."\n"; 
	while (&wantLineC($fh)) {
		my @ta=&splitL("\t", $_); 
		defined $back{$ta[0]} and die "repeat: $_\n"; 
		$back{$ta[0]} = [$cum, $ta[1]]; 
		defined $opts{'out_cumLen'} and print {$glob{'fh_oCumLen'}} join("\t", $ta[0], $ta[1], $back{$ta[0]}[0]+1, $back{$ta[0]}[0]+$ta[1])."\n"; 
		$cum += ($ta[1] + $opts{'delt_len'}); 
	}
	close($fh); 
	defined $opts{'out_cumLen'} and close($glob{'fh_oCumLen'}); 
	return(\%back); 
}

