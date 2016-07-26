#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 

!@ARGV and die "perl $0   marker_loc2GM_man   apple.snp\n"; 

my $fn_gm  = shift; 
my $fn_snp = shift; 

my $fh_gm  = &openFH($fn_gm, '<'); 
my %gm_anchors; 
my %rc_anchors; 
while (&wantLineC($fh_gm)) {
	my @ta = &splitL("\t", $_); 
	$ta[0] eq 'MarkerID' and next; 
	push(@{$gm_anchors{$ta[1]}}, [ $ta[2], $ta[8] ]); # [ phyP, gmP ]
}
close($fh_gm); 

# Reformat anchors. 
for my $cid (keys %gm_anchors) {
	my $min_phyP = $gm_anchors{$cid}[0][0]; 
	my $min_gmP  = $gm_anchors{$cid}[0][1]; 
	my $max_phyP = $gm_anchors{$cid}[-1][0]; 
	my $max_gmP  = $gm_anchors{$cid}[-1][1]; 
	my $d1 = ( $max_phyP > $min_phyP ) ? 1 : -1 ; 
	my $d2 = ( $max_gmP  > $min_gmP  ) ? 1 : -1 ; 
	my $strand = $d1 * $d2; 

	($min_phyP > $max_phyP) and ($min_phyP, $max_phyP) = ($max_phyP, $min_phyP); 
	($min_gmP  > $max_gmP ) and ($min_gmP , $max_gmP ) = ($max_gmP , $min_gmP ); 

	my ($gmD_per_bp_H , $gmD_per_bp_T); 
	@{$gm_anchors{$cid}} = sort { $a->[0] <=> $b->[0] } @{$gm_anchors{$cid}}; 
	for ( my $i=1; $i< @{$gm_anchors{$cid}}; $i++ ) {
		# $gm_anchors{$cid}[0][0] == 1 and last; 
		my $dist = abs( $gm_anchors{$cid}[$i][1] - $gm_anchors{$cid}[0][1] ) ; 
		if ( $dist > 1 ) {
			$gmD_per_bp_H = $dist / ( $gm_anchors{$cid}[$i][0] - $gm_anchors{$cid}[0][0] ); 
			# unshift( @{$gm_anchors{$cid}}, [ 1, $gm_anchors{$cid}[0][1] - $strand * $gmD_per_bpH * ( $gm_anchors{$cid}[0][0] -1 ) ] ); 
			last; 
		}
	}

	for ( my $i=$#{$gm_anchors{$cid}}-1; $i >= 0; $i-- ) {
		my $dist = abs( $gm_anchors{$cid}[-1][1] - $gm_anchors{$cid}[$i][1] ); 
		if ( $dist > 1 ) {
			$gmD_per_bp_T = $dist / ( $gm_anchors{$cid}[-1][0] - $gm_anchors{$cid}[$i][0] ); 
			last; 
		}
	}

	my ($new_min_gmP, $new_max_gmP); 

	for (my $i=0; $i<@{$gm_anchors{$cid}}; $i++) {
		my ($gmS, $gmE); 
		my ($min, $max); 
		if ( $i == 0 and $gm_anchors{$cid}[$i][0] > 1 ) {
			$gmS = $gm_anchors{$cid}[$i][1] - $strand * $gmD_per_bp_H * ($gm_anchors{$cid}[$i][0]-1); 
			$gmE = $gm_anchors{$cid}[$i][1] ; 
			push( @{$rc_anchors{$cid}}, [ 1, $gm_anchors{$cid}[$i][0]-1, $gmS , $gmE ] ); # [ phyS, phyE, gmS, gmE ]
			($min, $max) = &mathSunhh::minmax( [$gmS, $gmE] ); 
			$new_min_gmP //= $min; $new_min_gmP > $min and $new_min_gmP = $min; 
			$new_max_gmP //= $max; $new_max_gmP < $max and $new_max_gmP = $max; 
		}
		if ( $i < $#{$gm_anchors{$cid}} ) {
			$gmS = $gm_anchors{$cid}[$i][1]; 
			$gmE = $gm_anchors{$cid}[$i+1][1]; 
			push( @{$rc_anchors{$cid}}, [ $gm_anchors{$cid}[$i][0], $gm_anchors{$cid}[$i+1][0], $gmS , $gmE ] ); 
			($min, $max) = &mathSunhh::minmax( [$gmS, $gmE] ); 
			$new_min_gmP //= $min; $new_min_gmP > $min and $new_min_gmP = $min; 
			$new_max_gmP //= $max; $new_max_gmP < $max and $new_max_gmP = $max; 
		}
		if ( $i == $#{$gm_anchors{$cid}} ) {
			$gmS = $gm_anchors{$cid}[$i][1]; 
			$gmE = $gm_anchors{$cid}[$i][1] + $strand * $gmD_per_bp_T * 100e6; 
			push( @{$rc_anchors{$cid}}, [ $gm_anchors{$cid}[$i][0], $gm_anchors{$cid}[$i][0] + 100e6 , $gmS , $gmE ] ); 
			($min, $max) = &mathSunhh::minmax( [$gmS, $gmE] ); 
			$new_min_gmP //= $min; $new_min_gmP > $min and $new_min_gmP = $min; 
			$new_max_gmP //= $max; $new_max_gmP < $max and $new_max_gmP = $max; 
		}
	}

	if ( $strand == 1 ) {
		for my $t1 ( @{$rc_anchors{$cid}} ) {
			$t1->[2] = $t1->[2] - $new_min_gmP ; 
			$t1->[3] = $t1->[3] - $new_min_gmP ; 
		}
	} elsif ( $strand == -1 ) {
		for my $t1 ( @{$rc_anchors{$cid}} ) {
			$t1->[2] = $new_max_gmP - $t1->[2]; 
			$t1->[3] = $new_max_gmP - $t1->[3]; 
		}
	} else {
		die "strand=|$strand|\n"; 
	}

	&tsmsg("$cid\t$new_min_gmP\t$new_max_gmP\t$strand\n"); 
}

# Set Genetic map position for SNP table. 
my $fh_snp = &openFH( $fn_snp , '<' ); 
print STDOUT join("\t", qw/chr pos cM/)."\n"; 
my %cnt = ( 'cntN_base' => 0 , 'cntN_step' => 1e6 ); 
while (&wantLineC($fh_snp)) {
	&fileSunhh::log_section( $. , \%cnt ) and &tsmsg("[Msg] $. line.\n"); 
	m/^([^\t]+)\t([^\t]+)/ or die "$_\n"; 
	my ($chrID, $chrP) = ($1, $2); 
	$chrID eq 'chr' and next; 
	$chrID eq 'plast' and next; 
	$chrID eq 'mito' and next; 
	defined $rc_anchors{$chrID} or die "No data for [$chrID]\n"; 
	my $gmP; 
	for my $t1 ( @{$rc_anchors{$chrID}} ) {
		$chrP < $t1->[0] and next; 
		$chrP > $t1->[1] and next; 
		$gmP = $t1->[2] + ($t1->[3]-$t1->[2]) * ($chrP - $t1->[0]) / ( $t1->[1]-$t1->[0] ); 
		last; 
	}
	defined $gmP or die "failed for [$chrID , $chrP]\n"; 
	print STDOUT join("\t", $chrID, $chrP, $gmP)."\n"; 
}
close($fh_snp); 

&tsmsg("[Rec] Finish [$0]\n"); 


