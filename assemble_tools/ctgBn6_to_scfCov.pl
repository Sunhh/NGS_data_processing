#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 
use fileSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"minIdent:f", # 0.9 
	"log_lineN:i", # -1 
); 

$opts{'minIdent'} //= 0.9; 

my $help_txt = <<HH; 

perl $0 Gourd969_pilon_V1.ctg2scf.agp  blastn_100k_ctg.txt

-help 

-minIdent     [0.9] 0-1

-log_lineN    [-1]

HH

!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt);
$opts{'log_lineN'} //= -1; 

# [Sunhh@Falcon temp]$ head -3 blastn_100k_ctg.txt
# scaffold36_pilon_1      scaffold36_pilon_1      100.00  10712   0       0       1       10712   1       10712   0.0     19318   10712   10712   plus
# scaffold36_pilon_1      scaffold1_pilon_10      98.52   7072    93      5       1527    8588    23436   16367   0.0     12269   10712   63189   minus
# scaffold36_pilon_1      scaffold1_pilon_10      96.88   1027    31      1       6216    7242    46662   45637   0.0      1705   10712   63189   minus
#
# [Sunhh@wwz chk_cds_redundance]$ head -4 P1Genom_V1p2.ctg2scf.agp
# Cma_Scf00001    1       3213    1       W       Cma_Scf00001_1  1       3213    +
# Cma_Scf00001    3214    3396    2       N       183     scaffold        yes     paired-ends
# Cma_Scf00001    3397    5960    3       W       Cma_Scf00001_2  1       2564    +
# Cma_Scf00001    5961    6037    4       N       77      scaffold        yes     paired-ends
#

my %info; 
my $f1 = shift; 
&tsmsg("[Msg] Loading AGP file [$f1]\n"); 
my %ctg2scf = %{ &fileSunhh::load_agpFile( $f1 ) }; 
for my $ctgID (%ctg2scf) {
	for my $tb ( @{$ctg2scf{$ctgID}} ) {
		my @tc = @$tb; 
		$info{$tc[2]}{'noN'} += ($tc[1]-$tc[0]+1); 
	}
}
&tsmsg("[Msg] AGP file loaded.\n"); 

my %cnt; 
$cnt{'cntN_step'} = $opts{'log_lineN'}; 
my %blk; 
while (<>) {
	&fileSunhh::log_section( $. , \%cnt ) and &tsmsg("[Msg] Processing $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] >= 100 * $opts{'minIdent'} or next; 
	defined $ctg2scf{$ta[0]} or &stopErr("[Err] No contig_ID [$ta[0]] found.\n"); 
	my @r1s = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID'=>$ta[0], 'qryPos'=>$ta[6] ); 
	my @r1e = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID'=>$ta[0], 'qryPos'=>$ta[7] ); 
	$r1s[0][0] eq $r1e[0][0] or &stopErr("[Err] Bad scaffold found [$r1s[0][0]] and [$r1e[0][0]] : $_\n"); 
	my @r2s; 
	if (defined $ctg2scf{$ta[1]}) {
		@r2s = $ms_obj->switch_position( 'qry2ref'=>\%ctg2scf , 'qryID'=>$ta[1], 'qryPos'=>$ta[8] ); 
	} else {
		# defined $ctg2scf{$ta[1]} or &stopErr("[Err] No contig_ID [$ta[1]] found.\n"); 
		@r2s = ([$ta[1], $ta[8], '+']); 
	}
	my $scfID1 = $r1s[0][0]; 
	my $scfID2 = $r2s[0][0]; 
	$scfID1 eq $scfID2 and next; 
	my @se1 = ( $r1s[0][1], $r1e[0][1] ); 
	$se1[0] > $se1[1] and @se1[0,1] = @se1[1,0]; 
	push(@{$blk{$scfID1}{$scfID2}}, [$se1[0], $se1[1]]); 
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


