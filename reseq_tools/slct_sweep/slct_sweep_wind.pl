#!/usr/bin/perl
# Select windows for given selecting sweep rules. 
#  1. Use quantile value to find basic windows with smaller/equal ratio; 
#  2. Slide several windows as a group and detect good group fitting the required rules: 
#    2.1 The group length (window#) is set as 5; 
#    2.2 There should be at least three windows selected as basic windows; 
#    2.3 The edge of final group should be basic windows. 
#    In this way, windows 'Good Good Good Bad Bad Good Good' will fall into a same group, 
#      'Good Good Bad Good Bad Bad Good Good Good' will fall into a same group too. 
#      'Good Good Bad Good Bad Bad Bad Good Good Good' will fall into two different groups. 
#  3. Output all valid groups, as well as windows not covered by any group, as selected regions. 
#  4. Not sure if I should extend the selected regions with one window. 

use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use Statistics::Descriptive; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"inRatioFile:s", # 'apple_noN_w10ks10k.p1.all_grp1_grp2a_grp2b_grp2_grp3_grp4_grp5_grp2p3.PIavg.compare'
	"qtCutoff:f", # 0.01, percentile in perl Statistics::Descriptive module; 
	"qtFromLow!", 
	"grpLen:i", # 5, number of windows in a group checked. 
	"grpGood:i", # 3, number of good (selected) windows in a group checked. 
	"grpExtend:i", # 0, number of windows to extend group. 
	"slct_colN:i", # 5 
	"bpCnt_colN:i", # 4
); 

$opts{'grpLen'}     //= 5; 
$opts{'grpGood'}    //= 3; 
$opts{'qtCutoff'}   //= 0.03; 
$opts{'grpExtend'}  //= 0; 
$opts{'slct_colN'}  //= 5; 
$opts{'bpCnt_colN'} //= 4; 

$opts{'help'} and usage(); 
defined $opts{'inRatioFile'} or usage(); 
sub usage {
	print <<HH; 
################################################################################
# perl $0 -inRatioFile apple_noN_w10ks10k.p1.PIavg.compareRatioCol
# 
# -help
# 
# -qtCutoff        [$opts{'qtCutoff'}] percentile in perl Statistics::Descriptive module; 
# -qtFromLow       [Boolean] Count lower quantile if given. 
# -grpLen          [5] number of windows in a group checked
# -grpGood         [3] number of good (selected) windows in a group checked. 
# -grpExtend       [0] number of windows to extend group
#
# -slct_colN       [5] 
# -bpCnt_colN      [4]
################################################################################
HH
	exit 1; 
}

my $inRatioFh = &fileSunhh::openFH( $opts{'inRatioFile'}, '<' ); 
# ChrID   WindS   WindE   WindL   BpCnt C9_to_C6_min0.01
# chr1    1       10000   10000   9728  0.616125111088503
my %windInfor; 
while (<$inRatioFh>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1 and $ta[0] =~ m!^(ChrID|ChromID)$!i) {
		next; 
	}
	push(@{$windInfor{$ta[0]}{'line'}}, [@ta]); 
	$opts{'bpCnt_colN'} >= 0 and $ta[ $opts{'bpCnt_colN'} ] > 0 and push(@{$windInfor{$ta[0]}{'ratio_array'}}, $ta[ $opts{'slct_colN'} ]); 
}

# Select basic good windows. 
my $stat = Statistics::Descriptive::Full->new();
for my $chrID (sort keys %windInfor) {
	$stat->add_data( @{$windInfor{$chrID}{'ratio_array'}} ); 
	# &tsmsg("[Msg] @{$windInfor{$chrID}{'ratio_array'}}[0,1,2,3]\n"); 
}
if ( $opts{'qtFromLow'} ) {
	my $perc_tile = $stat->percentile( 100 * $opts{'qtCutoff'} ); 
	&tsmsg("[Rec] Percentile for $opts{'qtCutoff'} (from lower) is $perc_tile\n"); 
	for my $chrID ( sort keys %windInfor ) {
		for (my $i=0; $i<@{ $windInfor{$chrID}{'line'} }; $i++) {
			$opts{'bpCnt_colN'} >= 0 and do { $windInfor{$chrID}{'line'}[$i][ $opts{'bpCnt_colN'} ] > 0 or next }; 
			$windInfor{$chrID}{'line'}[$i][ $opts{'bpCnt_colN'} ] <= $perc_tile and push(@{$windInfor{$chrID}{'goodIdx'}}, $i); 
		}
	}
} else {
	my $perc_tile = $stat->percentile( 100 - 100 * $opts{'qtCutoff'} ); 
	&tsmsg("[Rec] Percentile for $opts{'qtCutoff'} (from higher) is $perc_tile\n"); 
	for my $chrID ( sort keys %windInfor ) {
		for (my $i=0; $i<@{ $windInfor{$chrID}{'line'} }; $i++) {
			$opts{'bpCnt_colN'} >= 0 and do { $windInfor{$chrID}{'line'}[$i][ $opts{'bpCnt_colN'} ] > 0 or next; }; 
			$windInfor{$chrID}{'line'}[$i][ $opts{'bpCnt_colN'} ] =~ m/^[\d.]+$/i or next; 
			$windInfor{$chrID}{'line'}[$i][ $opts{'bpCnt_colN'} ] >= $perc_tile and push(@{$windInfor{$chrID}{'goodIdx'}}, $i); 
		}
	}
}

# Search for good groups. 
for my $chrID (sort keys %windInfor) {
	defined $windInfor{$chrID}{'goodIdx'} or next; 
	&group_good_winds( $windInfor{$chrID}, $opts{'grpLen'}, $opts{'grpGood'}, $opts{'grpExtend'} ); 
}

# Output good groups and windows. 
for my $chrID (sort keys %windInfor) {
	&output_winds( $windInfor{$chrID} ); 
}

################################################################################
# Sub-routines. 
################################################################################
sub output_winds {
	my ($wind_hr) = @_; 
	$wind_hr->{'grp_se'} //= []; 
	my @out_se_i = @{$wind_hr->{'grp_se'}}; 
	$wind_hr->{'goodIdx'} //= []; 
	for (my $i=0; $i<@{$wind_hr->{'goodIdx'}}; $i++) {
		defined $wind_hr->{'inGrp'}{ $wind_hr->{'goodIdx'}[$i] } and next; 
		push(@out_se_i, [$wind_hr->{'goodIdx'}[$i], $wind_hr->{'goodIdx'}[$i]]); 
	}
	for my $ar ( sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @out_se_i ) {
		my @loc1 = @{$wind_hr->{'line'}[$ar->[0]]} ; 
		my @loc2 = @{$wind_hr->{'line'}[$ar->[1]]} ; 
		my $sum_bpCnt = 0; 
		for my $idx ( $ar->[0] .. $ar->[1] ) {
			if ( $opts{'bpCnt_colN'} >= 0 ) {
				$sum_bpCnt += $wind_hr->{'line'}[ $idx ][ $opts{'bpCnt_colN'} ]; 
			} else {
				$sum_bpCnt += ( $wind_hr->{'line'}[ $idx ][2] - $wind_hr->{'line'}[ $idx ][1] + 1 ); 
			}
		}
		print STDOUT join("\t", $loc1[0], $loc1[1], $loc2[2], $loc2[2]-$loc1[1]+1, $sum_bpCnt)."\n"; 
	}

	return ; 
}# output_winds() 

sub group_good_winds {
	my ($hr, $len, $goodN, $sideEnlarge) = @_; 
	$len //= 5; 
	$goodN //= 3; 
	$sideEnlarge //= 0; 
	my %h = %$hr; 
	my @ti = @{$h{'goodIdx'}}; 
	my @grp; 
	for (my $i=0; $i<@ti; $i++) {
		my $good_in_grp = 1; 
		my $startIdx = $ti[$i]; 
		my $endIdx = $startIdx + $len - 1; 
		for (my $j=$i+1; $j<@ti and $ti[$j] <= $endIdx; $j++) {
			$good_in_grp ++; 
		}
		if ( $good_in_grp >= $goodN ) {
			# This is a good group. 
			&extend_grp(\@grp, $startIdx, $endIdx); 
		}
	}
	&chop_grp(\@grp, \@ti); # Chop groups to limit edge. 
	&enlarge_grp(\@grp, $sideEnlarge); 
	
	@{$hr->{'grp_se'}} = @grp; 
	for my $tr1 (@grp) {
		my ( $s, $e ) = @$tr1; 
		for my $idx ( $s .. $e ) {
			$hr->{'inGrp'}{$idx} = 1; 
		}
	}
	return ; 
}# group_good_winds

sub extend_grp {
	my ($ar, $s, $e) = @_; 
	if (@$ar > 0) {
		if ( $ar->[-1][1] >= $s ) {
			$ar->[-1][1] < $e and $ar->[-1][1] = $e; 
		} else {
			push(@$ar, [$s, $e]); 
		}
	} else {
		push(@$ar, [$s, $e]); 
	}
}# sub extend_grp () 

sub chop_grp {
	my ($grp_ar, $ti_ar) = @_; 
	my %good; 
	for (@$ti_ar) {
		$good{$_} = 1; 
	}
	for my $ar (@$grp_ar) {
		my ($s, $e) = @$ar; 
		my ($new_s, $new_e); 
		for my $idx ($s .. $e) {
			defined $good{$idx} or next; 
			defined $new_s or $new_s = $idx; 
			$new_e = $idx; 
		}
		@$ar = ($new_s, $new_e); 
	}
	return; 
}# chop_grp ()

sub enlarge_grp {
	my ($grp_ar, $enlargeN) = @_; 
	for my $ar (@$grp_ar) {
		my $new_s = $ar->[0] - $enlargeN; 
		my $new_e = $ar->[1] + $enlargeN; 
		$new_s >= 0 or $new_s = 0; 
		@$ar = ($new_s, $new_e); 
	}
	return ; 
}# enlarge_grp() 

#$opts{'grpLen'} //= 5; 
#$opts{'grpGood'} //= 3; 


