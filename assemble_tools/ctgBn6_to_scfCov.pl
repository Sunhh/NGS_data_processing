#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"minIdent:f", # 0.9 
); 

$opts{'minIdent'} //= 0.9; 

my $help_txt = <<HH; 

perl $0 Gourd969_pilon_V1.ctg.fa.kl  blastn_100k_ctg.txt

-help 

-minIdent     [0.9] 0-1

HH

!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt);


# [Sunhh@Falcon temp]$ head -3 blastn_100k_ctg.txt
# scaffold36_pilon_1      scaffold36_pilon_1      100.00  10712   0       0       1       10712   1       10712   0.0     19318   10712   10712   plus
# scaffold36_pilon_1      scaffold1_pilon_10      98.52   7072    93      5       1527    8588    23436   16367   0.0     12269   10712   63189   minus
# scaffold36_pilon_1      scaffold1_pilon_10      96.88   1027    31      1       6216    7242    46662   45637   0.0      1705   10712   63189   minus
#
# [Sunhh@Falcon temp]$ head Gourd969_pilon_V1.ctg.fa.kl
# key     len
# scaffold36_pilon_1      10712
# scaffold36_pilon_2      21776
# scaffold36_pilon_3      5401
#

my (%info); 
my $f1 = shift; 
open F1,'<',"$f1" or die; 
while (<F1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] =~ m/^key$/i and next; 
	my ($scfID) = &sep_ctgID($ta[0]); 
	$info{$scfID}{'noN'} += $ta[1]; 
}
close F1; 

my %blk; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] >= 100 * $opts{'minIdent'} or next; 
	my ($scfID1) = &sep_ctgID( $ta[0] ); 
	my ($scfID2) = &sep_ctgID( $ta[1] ); 
	$scfID1 eq $scfID2 and next; 
	$ta[6] > $ta[7] and @ta[6,7] = @ta[7,6]; 
	push(@{$blk{$scfID1}{$scfID2}}, [$ta[6], $ta[7]]); 
}

print STDOUT join("\t", qw/ScfID LenInCtg CovLenInScf CovPercInScf CovByScf CovLenSpan CovPercSpan/)."\n"; 
for my $tk1 (keys %blk) {
	my %maxSum; 
	my @tb; 
	for my $tk2 (keys %{$blk{$tk1}}) {
		$blk{$tk1}{$tk2} = $ms_obj->mergeLocBlk( $blk{$tk1}{$tk2}, 'dist2join'=>1 ); 
		push(@tb, $blk{$tk1}{$tk2}); 
		my $tsum = 0; 
		for my $ar1 ( @{ $blk{$tk1}{$tk2} } ) {
			$tsum += $ar1->[1]-$ar1->[0]+1; 
		}
		$maxSum{'tk2'} //= $tk2; 
		$maxSum{'sum'} //= $tsum; 
		if ($maxSum{'sum'} < $tsum) {
			$maxSum{'tk2'} = $tk2; 
			$maxSum{'sum'} = $tsum; 
		}
	}

	my $tc = $ms_obj->mergeLocBlk( @tb, 'dist2join'=>1 ); 
	my $sum_span = 0; 
	for my $ar2 (@$tc) {
		$sum_span += $ar2->[1]-$ar2->[0]+1; 
	}

	print STDOUT join("\t", 
	  $tk1, 
	  $info{$tk1}{'noN'}, 
	  $maxSum{'sum'}, 
	  sprintf("%0.2f", $maxSum{'sum'}/$info{$tk1}{'noN'}*100), 
	  $maxSum{'tk2'}, 
	  $sum_span, 
	  sprintf("%0.2f", $sum_span / $info{$tk1}{'noN'}*100)
	)."\n"; 
}

sub sep_ctgID {
	$_[0] =~ m/^(\S+)_(\d+)$/ or die "bad ID [$_[0]]\n"; 
	return($1, $2); 
}


