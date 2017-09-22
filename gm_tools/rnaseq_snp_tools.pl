#!/usr/bin/perl -w 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
use SNP_tbl; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	# Tasks : 
	"task:s", 
	
	# task: cnt_pop_allele() 
	"vcf_tab:s", # Input file coming from "vcf_to_tab"; The first three columns are #CHROM , POS and REF . 
		"P1_colN:i", "P2_colN:i", # column number for two parents (P1 and P2) 
		"offs_colN:s", 
		"infer_PG!", 
		"keep_abnormal!", 
	
	# task: get_abh() 
	#  Same parameters to task:cnt_pop_allele, 
	#  but the following don't work: 'keep_abnormal'
	
	# task: cnt_abh() 
	
	# task: blk_abh() 
	"abh_tab:s",    # Output of task:get_abh; 
	"len_jnBlk:i",  # 1000000; Length within which adjacent blocks with same inferred genotype are joined. 
	"skipU_jnBlk!", # Boolean; If given, genotype 'u' will be ignored when joining blocks. 
	"window_geno!", # This should be executed before -smooth_geno. 
		"len_window:i",     # 200e3; 
		"sN2_window:i",     # 10; At least this number of siteN2 in window. 
		"window_round:i",   # 10; Maximum number of rounds to window genotypes.
		"window_hRatio:f",  # 0.2; If having this ratio of heterozygous sites, this window is heterozygous; 
		"window_aRatio:f",  # 0.3; If having this ratio of 'a' sites, this window has 'a' allele; 
		"window_bRatio:f",  # 0.3; If having this ratio of 'b' sites, this window has 'b' allele; 
	"smooth_geno!", 
		"smooth_round:i",   # 10; 
		"len_noX:i",        # 1000 ; Length within which there should be no crossing; 
		"len_singleP:i",    # 200  ; Length within which we count at most one SNP site (marker); 
		"len_dbX:i",        # 100e3; Length within which we check double crossing problem; 
		"fixGeno_siteN:i",  # 10   ; If there are so many sites supporting this genotype, we don't change it at all (Except genotype 'u'); 
		"bigger_copyN:f", "bigger_miniN:i", "bigger_diffN:i",   
		# 1.5, 3, 3  ; siteN_1 is bigger than siteN_2 if all the followings three are true : 
		#    siteN_1 >= $bigger_copyN * siteN_2 
		#    siteN_1 >= $bigger_miniN 
		#    siteN_1 >= $bigger_diffN + siteN_2 
		# Or else, siteN_1 is similar to siteN_2; 
	
	# task: rejoin_abhblk() 
	"abhblk_tab:s",  # Output of task:blk_abh; 
	"min_isiteN2:i", # 0 
	# "skipU_jnBlk!"  still works. 
	# "len_jnBlk:i"   still works. 
	# "len_singleP:i" still works. 
	
); 

my %glob; 
&initial_glob(); 
&mathSunhh::_addHash('toH' => \%glob, 'fromH' => \%opts, 'replaceExist' => 1); 

$glob{'task'} eq 'NA' and &usage("\n"); 
$glob{'help'} and &usage("\n"); 

if ( $glob{'task'} =~ m!^cnt_pop_allele$!i ) {
	&cnt_pop_allele(); 
} elsif ( $glob{'task'} =~ m!^get_abh$!i ) {
	&get_abh(); 
} elsif ( $glob{'task'} =~ m!^blk_abh$!i ) {
	&blk_abh(); 
} elsif ( $glob{'task'} =~ m!^rejoin_abhblk$!i ) {
	&rejoin_abhblk(); 
}


####################################################################################################
#       Functional sub-routines 
####################################################################################################

sub rejoin_abhblk {
	&chk_para('rejoin_abhblk'); 
	
	my $jnBlkH = &_load_abhblk( $glob{'abhblk_tab'} ); 
	
	if ( $glob{'min_isiteN2'} > 0 ) {
		for my $indvID ( sort { $jnBlkH->{$a}{'order'} <=> $jnBlkH->{$b}{'order'} } keys %$jnBlkH) {
			&tsmsg("[Msg] Trimming $jnBlkH->{$indvID}{'order'}-th [$indvID]\n"); 
			my @filtered_blk; 
			for my $tb ( @{$jnBlkH->{$indvID}{'blk'}} ) {
				if ( &_is_same_genotype( $tb->[5][0][1], 'u' ) ) {
					push( @filtered_blk, $tb ); 
				} elsif ( $tb->[5][0][2] >= $glob{'min_isiteN2'} ) {
					push( @filtered_blk, $tb ); 
				} else {
#			$tb->[0] eq 'SL2.40ch12' and die "@$tb\n"; 
					; 
				}
			}
			$jnBlkH->{$indvID}{'blk'} = \@filtered_blk; 
		}
	}
	
	for my $indvID (sort { $jnBlkH->{$a}{'order'} <=> $jnBlkH->{$b}{'order'} } keys %$jnBlkH) {
		&tsmsg("[Msg] Joining $jnBlkH->{$indvID}{'order'}-th [$indvID]\n"); 
		my $j2 = &_join_blk_i( $jnBlkH->{$indvID}, 0, $glob{'len_jnBlk'}, $glob{'skipU_jnBlk'} ); 
		&output_blkH_i( $j2 ); 
	}
	
	return; 
}# rejoin_abhblk() 

#	# task: blk_abh() 
#	"len_jnBlk:i", # 1000000; Length within which adjacent blocks with same inferred genotype are joined. 
#	"skipU_jnBlk!",# Boolean; If given, genotype 'u' will be ignored when joining blocks. 
#	"smooth_geno!", 
#		"len_singleP:i",    # 200  ; Length within which we count at most one SNP site (marker); 
#		"len_noX:i",        # 1000 ; Length within which there should be no crossing; 
#		"len_dbX:i",        # 100e3; Length within which we check double crossing problem; 
#		"fixGeno_siteN:i",  # 10   ; If there are so many sites supporting this genotype, we don't change it at all (Except genotype 'u'); 
#		"bigger_copyN:f", "bigger_miniN:i", "bigger_diffN:i",   
#		# 1.5, 3, 3  ; siteN_1 is bigger than siteN_2 if all the followings three are true : 
#		#    siteN_1 >= $bigger_copyN * siteN_2 
#		#    siteN_1 >= $bigger_miniN 
#		#    siteN_1 >= $bigger_diffN + siteN_2 
#		# Or else, siteN_1 is similar to siteN_2; 
sub blk_abh {
	&chk_para('blk_abh'); 
	my $blkH = &_load_abh( $glob{'abh_tab'} ); 
	
	unless ( $glob{'smooth_geno'} or $glob{'window_geno'} ) {
		# If we don't require to 'smooth_geno', I can simply output abh regions by each individual. 
		#   If 'skipU_jnBlk' is true, the inferred 'U' will be ignored in calculating. 
		#   If pos_2-pos_1 <= 'len_jnBlk' and they are the same genotype, I'll join them; 
		for (my $indvI=0; $indvI<@{$blkH->{'offs'}}; $indvI++) {
			my $blkH_jn = &_join_blk_i( $blkH, $indvI, $glob{'len_jnBlk'}, $glob{'skipU_jnBlk'} ); 
			&output_blkH_i( $blkH_jn ); 
		}
		return; 
	}
	
	if ( $glob{'window_geno'} ) {
		for (my $indvI=0; $indvI<@{$blkH->{'offs'}}; $indvI++) {
			my $is_changed = 1; 
			my $cnt_changed_times = 0; 
			while ($is_changed) {
				$cnt_changed_times ++; 
				$cnt_changed_times > $glob{'window_round'} and last; 
				$is_changed = 0; 
				&tsmsg("[Msg] window_geno for $indvI-th [$blkH->{'offs'}[$indvI]] round [$cnt_changed_times].\n"); 
				BLK_J: 
				for (my $blkJ=0; $blkJ+1<@{$blkH->{'blk'}}; $blkJ++) {
					my $abh_curr  = $blkH->{'blk'}[$blkJ][5][$indvI][1]; 
					&_is_same_genotype( $abh_curr , 'u' ) and next BLK_J; 
					my $chID_curr = $blkH->{'blk'}[$blkJ][0]; 
					my $blkJ_next = $blkJ+1; 
					$blkH->{'blk'}[$blkJ_next][0] eq $chID_curr or next BLK_J; 
					BLK_J_Next: 
					for ( ; $blkJ_next < @{$blkH->{'blk'}}; $blkJ_next++ ) {
						if ( $blkH->{'blk'}[$blkJ_next][0] ne $chID_curr ) {
							$blkJ = $blkJ_next-1; 
							next BLK_J; 
						}
						$blkH->{'blk'}[$blkJ_next][1] - $blkH->{'blk'}[$blkJ][2] <= $glob{'len_window'} or next BLK_J; 
						&_is_same_genotype( $blkH->{'blk'}[$blkJ_next][5][$indvI][1], 'u' ) and next BLK_J_Next; 
						&_is_same_genotype( $abh_curr, $blkH->{'blk'}[$blkJ_next][5][$indvI][1] ) and do { $blkJ = $blkJ_next-1; next BLK_J; }; 
						last BLK_J_Next; 
					}# BLK_J_Next : for 
					$blkJ_next < @{$blkH->{'blk'}} or last BLK_J; # No difference found in all following sites. 
					my $is_changed_t = &window_blkJ_it( $blkH, $indvI, $blkJ, 0, $blkJ_next ); 
					if ( $is_changed_t == 1 ) {
						$is_changed = 1; 
						for (; $blkJ>=0; $blkJ--) {
							&_is_same_genotype($blkH->{'blk'}[$blkJ][5][$indvI][1], 'u') and next; 
							$blkH->{'blk'}[$blkJ_next][0] eq $blkH->{'blk'}[$blkJ][0] or last; 
							$blkH->{'blk'}[$blkJ_next][1] - $blkH->{'blk'}[$blkJ][2] <= $glob{'len_window'} or last ;
						}
						$blkJ >= 0 or $blkJ = 0; 
						$blkJ--; 
						next BLK_J; 
					}
				}# for:BLK_J; for (my $blkJ=0; $blkJ+1<@{$blkH->{'blk'}}; $blkJ++) 
			}# while ($is_changed) 
		}
	}# if ( $glob{'window_geno'} )
	
	# Now let's try to 'smooth_geno' for each offspring individual!!! 
	if ( $glob{'smooth_geno'} ) {
		for (my $indvI=0; $indvI<@{$blkH->{'offs'}}; $indvI++) {
			# Smooth from each block 
			my $is_changed = 1; 
			my $cnt_changed_times = 0; 
			while ($is_changed) {
				$cnt_changed_times++; 
				$cnt_changed_times > $glob{'smooth_round'} and last; 
				$is_changed = 0; 
				&tsmsg("[Msg] Processing $indvI-th [$blkH->{'offs'}[$indvI]] round [$cnt_changed_times].\n"); 
				for (my $blkJ=0; $blkJ+1<@{$blkH->{'blk'}}; $blkJ++) {
					my $abh_curr  = $blkH->{'blk'}[$blkJ][5][$indvI][1]; 
					my $chID_curr = $blkH->{'blk'}[$blkJ][0]; 
					&_is_same_genotype( $abh_curr, 'u' ) and next; 
					my $blkJ_next = $blkJ+1; 
					my $abh_next; 
					for (; $blkJ_next < @{$blkH->{'blk'}}; $blkJ_next++) {
						&_is_same_genotype( $blkH->{'blk'}[$blkJ_next][5][$indvI][1], 'u') and next; 
						if ( $chID_curr eq $blkH->{'blk'}[$blkJ_next][0] ) {
							if ( &_is_same_genotype( $abh_curr, $blkH->{'blk'}[$blkJ_next][5][$indvI][1] ) ) {
								$blkJ = $blkJ_next; 
							} else {
								$abh_next  = $blkH->{'blk'}[$blkJ_next][5][$indvI][1]; 
								last; 
							}
						} else {
							last; 
						}
					}
					defined $abh_next        or  do { $blkJ = $blkJ_next-1; next; }; 
					# &tsmsg("[Msg] Processing $indvI-th $blkH->{'offs'}[$indvI] at [$blkJ: @{$blkH->{'blk'}[$blkJ]}[0,1,2] $abh_curr] compared to [$blkJ_next: @{$blkH->{'blk'}[$blkJ_next]}[0,1,2] $abh_next]\n"); 
					my $is_changed_t = &smooth_blkJ_it( $blkH, $indvI, $blkJ, 0 ); # input format : ( &_load_abh(), $offspring_indexI, $block_indexJ, $changed_times_atJ )
					$is_changed_t == 1 and $is_changed = 1; 
				}
			}# while ($is_changed) 
		}
	}# if ( $glob{'smooth_geno'} ) 
	
	# Output smoothed blocks; 
	for (my $indvI=0; $indvI<@{$blkH->{'offs'}}; $indvI++) {
		my $blkH_jn = &_join_blk_i( $blkH, $indvI, $glob{'len_jnBlk'}, $glob{'skipU_jnBlk'} ); 
		&output_blkH_i( $blkH_jn ); 
	}
	
	return; 
}# blk_abh() 

sub window_blkJ_it{
	my ($blkH, $indvI, $blkJ, $runs_cnt, $blkJ_next) = @_; 
	$blkJ_next //= $blkJ+1; 
	my $is_changed = 1; 
	my $curr_cnt = $runs_cnt; 
	($is_changed, $curr_cnt) = &window_blkJ($blkH, $indvI, $blkJ, $curr_cnt, $blkJ_next); 
	my $has_been_changed = $is_changed;  
	while ( $is_changed ) {
		($is_changed, $curr_cnt) = &window_blkJ($blkH, $indvI, $blkJ, $curr_cnt, $blkJ_next); 
	}
	return($has_been_changed); 
}# window_blkJ_it() 

sub smooth_blkJ_it {
	my ($blkH, $indvI, $blkJ, $runs_cnt) = @_; 
	my $is_changed = 1; 
	my $curr_cnt = $runs_cnt; 
	($is_changed, $curr_cnt) = &smooth_blkJ($blkH, $indvI, $blkJ, $curr_cnt); 
	my $has_been_changed = $is_changed; 
	while ( $is_changed ) {
		($is_changed, $curr_cnt) = &smooth_blkJ($blkH, $indvI, $blkJ, $curr_cnt); 
	}
	return ($has_been_changed); 
}# smooth_blkJ_it() 

# ($is_changed, $curr_cnt) = &window_blkJ($blkH, $indvI, $blkJ, $curr_cnt, $blkJ_next); 
sub window_blkJ {
	my ($blkH, $indvI, $blkJ_0, $runs_cnt, $blkJ_1) = @_; 
	$runs_cnt //= 0; 
	$runs_cnt ++; 
	if ($runs_cnt > 100) {
		&tsmsg("[Wrn] Too many times [$runs_cnt-1] for window [offspring_indexI=$indvI] [$blkH->{'offs'}[$indvI]] [block_indexJ=$blkJ_0]\n"); 
		return(0, $runs_cnt); 
	}
	$blkJ_1 //= $blkJ_0 + 1; 
	my $abh_0 = $blkH->{'blk'}[$blkJ_0][5][$indvI][1]; 
	my $abh_1 = $blkH->{'blk'}[$blkJ_1][5][$indvI][1]; 
	&_is_same_genotype($abh_0, 'u') and return(0, $runs_cnt); 
	&_is_same_genotype($abh_1, 'u') and return(0, $runs_cnt); 
	&_is_same_genotype($abh_0, $abh_1) and return(0, $runs_cnt); 
	my $lwind = &type_wind($blkH, $indvI, $blkJ_0, 'left'); 
	my $rwind = &type_wind($blkH, $indvI, $blkJ_1, 'right'); 
	my $lwind_t = $lwind->{'wind_type'}; 
	my $rwind_t = $rwind->{'wind_type'}; 
	
	if ( &_is_same_genotype($lwind_t, 'u') ) {
		if ( &_is_same_genotype($rwind_t, 'u') ) {
			return( &replace_abh( $blkH, $indvI, [$blkJ_0, $blkJ_1], 'u', $glob{'fixGeno_siteN'}-1 ), $runs_cnt ); 
		} else {
			my $tc = 0; 
			unless ( &_is_same_genotype( $abh_0 , $rwind_t ) ) {
				&replace_abh( $blkH, $indvI, [$blkJ_0], 'u',                   $glob{'fixGeno_siteN'}-1 ) and $tc = 1; 
			}
			&replace_abh( $blkH, $indvI, [$blkJ_1], $rwind_t, $glob{'fixGeno_siteN'}-1 ) and $tc = 1; 
			return( $tc, $runs_cnt ); 
		}
	}
	if ( &_is_same_genotype($rwind_t, 'u') ) {
		return( &replace_abh( $blkH, $indvI, [$blkJ_1], $rwind_t, $glob{'fixGeno_siteN'}-1 ), $runs_cnt ); 
	}

#if ($blkH->{'offs'}[$indvI] eq 'F15-4231-rep2' and $blkH->{'blk'}[$blkJ_0][2] >= 1862380 and $blkH->{'blk'}[$blkJ_0][1] <= 1862380) {
if ($blkH->{'offs'}[$indvI] eq 'F15-3882-rep1' and $blkH->{'blk'}[$blkJ_0][0] eq 'SL2.40ch05' and $blkH->{'blk'}[$blkJ_0][2] == 1338805) {
	warn "abh_0=$abh_0;abh_1=$abh_1;\n"; 
	warn "lw=$lwind_t;rw=$rwind_t\n"; 
	warn "lw_p=$lwind->{'wind_abh'}[0][1] $lwind->{'wind_abh'}[-1][1]; lw_a=$lwind->{'wind_abhR'}{'a_r'} lw_b=$lwind->{'wind_abhR'}{'b_r'} lw_h=$lwind->{'wind_abhR'}{'h_r'}\n"; 
	warn "rw_p=$rwind->{'wind_abh'}[0][1] $rwind->{'wind_abh'}[-1][1]; rw_a=$rwind->{'wind_abhR'}{'a_r'} rw_b=$rwind->{'wind_abhR'}{'b_r'} rw_h=$rwind->{'wind_abhR'}{'h_r'}\n"; 
	warn "[@{$blkH->{'blk'}[$blkJ_0]}[0..3]]\n"; 
	die "\n"; 
}
	unless ( &_is_same_genotype($abh_0, $lwind_t) ) {
for (my $i=1000; $i>=1; $i--) {
	my $blkJ_new = $blkJ_0 - $i; 
	&tsmsg("[Err] Prev: @{$blkH->{'blk'}[$blkJ_new]}[0..3] @{$blkH->{'blk'}[$blkJ_new][5][$indvI]}\n"); 
}
		&stopErr(
			"[Err] Bad site[$abh_0 $abh_1] VS. wind[$lwind_t $rwind_t] lw[$lwind->{'wind_abh'}[0][1]..$lwind->{'wind_abh'}[-1][1]] rw[$rwind->{'wind_abh'}[0][1]..$rwind->{'wind_abh'}[-1][1]] in $blkH->{'offs'}[$indvI] [@{$blkH->{'blk'}[$blkJ_0]}[0..4] [@{$blkH->{'blk'}[$blkJ_0][5][$indvI]}]] @{$blkH->{'blk'}[$lwind->{'wind_E_i'}]}[0..3]\n"
		); 
	}
	if      ( !(&_is_same_genotype($abh_0, 'h')) and &_is_same_genotype($abh_1, 'h')    ) {
		if      ( &_is_same_genotype($rwind_t, 'h') ) {
			return( &replace_abh( $blkH, $indvI, [$rwind->{'wind_S_i'} .. $rwind->{'wind_E_i'}], $rwind_t, $glob{'fixGeno_siteN'}-1 ), $runs_cnt ); 
		} elsif ( &_is_same_genotype($rwind_t, $lwind_t) ) {
			return( &replace_abh( $blkH, $indvI, [$blkJ_1], $rwind_t, $glob{'fixGeno_siteN'}-1 ), $runs_cnt ) ; 
		} else {
			return( &replace_abh( $blkH, $indvI, [$blkJ_1], 'u',                   $glob{'fixGeno_siteN'}-1 ), $runs_cnt ) ; 
		}
	} elsif ( !(&_is_same_genotype($abh_0, 'h')) and !(&_is_same_genotype($abh_1, 'h'))  ) {
		if      ( &_is_same_genotype($rwind_t, 'h') ) {
			return( &replace_abh( $blkH, $indvI, [$rwind->{'wind_S_i'} .. $rwind->{'wind_E_i'}], $rwind_t, $glob{'fixGeno_siteN'}-1 ), $runs_cnt ); 
		} elsif ( &_is_same_genotype($rwind_t, $lwind_t) ) {
			return( &replace_abh( $blkH, $indvI, [$blkJ_1], $rwind_t, $glob{'fixGeno_siteN'}-1 ), $runs_cnt ) ; 
		} else {
			return( &replace_abh( $blkH, $indvI, [$rwind->{'wind_S_i'} .. $rwind->{'wind_E_i'}], $rwind_t, $glob{'fixGeno_siteN'}-1 ), $runs_cnt ) ; 
		}
	} elsif ( &_is_same_genotype($abh_0, 'h')    and !(&_is_same_genotype($abh_1, 'h')) ) {
		if ( !(&_is_same_genotype($rwind_t, 'h')) ) {
			if ( &_is_same_genotype($abh_1, $rwind_t) ) {
				### Edit here for s1-H s2-AB + w-H + w-AB 
				return( &replace_abh( $blkH, $indvI, [$rwind->{'wind_S_i'} .. $rwind->{'wind_E_i'}], $rwind_t, $glob{'fixGeno_siteN'}-1 ) , $runs_cnt ); 
			} else {
				return( &replace_abh( $blkH, $indvI, [$blkJ_1], 'u', $glob{'fixGeno_siteN'}-1 ), $runs_cnt ); 
			}
		} else {
			return( &replace_abh( $blkH, $indvI, [$rwind->{'wind_S_i'} .. $rwind->{'wind_E_i'}], 'h', $glob{'fixGeno_siteN'}-1 ), $runs_cnt ); 
		}
	}
	return(0, $runs_cnt); 
}# window_blkJ() 


sub smooth_blkJ {
	my ($blkH, $indvI, $blkJ_0, $runs_cnt) = @_; 
	$runs_cnt //= 0; 
	$runs_cnt ++; 
	if ($runs_cnt > 100) {
		# There might be some problem if this function is iterated too many times. 
		#   It might be a good idea to output the local section. 
		&tsmsg("[Wrn] Too many times [$runs_cnt-1] for smoothing [offspring_indexI=$indvI][$blkH->{'offs'}[$indvI]] [block_indexJ=$blkJ_0]\n"); 
		return(0, $runs_cnt); 
	}
	my $sP_0   = $blkH->{'blk'}[$blkJ_0][1]; 
	my $eP_0   = $blkH->{'blk'}[$blkJ_0][2]; 
	my $chID_1 = $blkH->{'blk'}[$blkJ_0][0]; 
	my $sP_1   = $blkH->{'blk'}[$blkJ_0][1]; 
	my $eP_1   = $blkH->{'blk'}[$blkJ_0][2]; 
	my $blkJ_1s= $blkJ_0; 
	my $blkJ_1e= $blkJ_0; 
	my ($n1,    $n2,    $n3,    $n4); 
	my ($abh_1, $abh_2, $abh_3, $abh_4); 
	my (@blkJ_idx_1, @blkJ_idx_2, @blkJ_idx_3, @blkJ_idx_4); 
	$abh_1 = $blkH->{'blk'}[$blkJ_0][5][$indvI][1]; 
	&_is_same_genotype($abh_1, 'u') and return(0, $runs_cnt); 
	
	# Extend to the left for n1; 
	for (my $i=$blkJ_0-1; $i>=0; $i--) {
		$chID_1 eq $blkH->{'blk'}[$i][0] or last; 
		$sP_1 > $blkH->{'blk'}[$i][2] or &stopErr("[Err] left error [$chID_1 $sP_1 back to $blkH->{'blk'}[$i][2]].\n"); 
		$sP_1 - $blkH->{'blk'}[$i][2] <= $glob{'len_dbX'} or last; 
		&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; # I am trying to include same blocks spanning 'u' here. 
		&_is_same_genotype($abh_1, $blkH->{'blk'}[$i][5][$indvI][1]) or last; 
		$sP_1   = $blkH->{'blk'}[$i][1]; 
		$blkJ_1s= $i; 
	}
	# Extend to the right for n1; 
	for (my $i=$blkJ_0+1; $i<@{$blkH->{'blk'}}; $i++) {
		$chID_1 eq $blkH->{'blk'}[$i][0] or last; 
		$eP_1 < $blkH->{'blk'}[$i][1] or &stopErr("[Err] right error n1 [$chID_1 $eP_1 forward to $blkH->{'blk'}[$i][1]]\n"); 
		$blkH->{'blk'}[$i][1] - $eP_1 <= $glob{'len_dbX'} or last; 
		&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
		&_is_same_genotype($abh_1, $blkH->{'blk'}[$i][5][$indvI][1]) or last; 
		$eP_1   = $blkH->{'blk'}[$i][2]; 
		$blkJ_1e= $i; 
	}
	$eP_1 - $eP_0 <= $glob{'len_dbX'} or return(0, $runs_cnt); 
	@blkJ_idx_1 = ($blkJ_1s .. $blkJ_1e); 
	$n1 = &get_blkIndvN2( $blkH, $indvI, \@blkJ_idx_1, $abh_1 ); 
	# $blkH->{'offs'}[$indvI] eq 'F15-4231-rep1' and warn "$chID_1 [ $blkH->{'blk'}[$blkJ_1s][1] - $blkH->{'blk'}[$blkJ_1s][2] ]: @blkJ_idx_1\n"; 
	
	# Extend to the right for n2 from n1; 
	my ($sP_2, $eP_2, $blkJ_2s, $blkJ_2e); 
	for (my $i=$blkJ_1e+1; $i<@{$blkH->{'blk'}}; $i++) {
		$chID_1 eq $blkH->{'blk'}[$i][0] or last; 
		if (defined $abh_2) {
			$blkH->{'blk'}[$i][1] - $eP_2 <= $glob{'len_dbX'} or last; 
			&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
			&_is_same_genotype($abh_2, $blkH->{'blk'}[$i][5][$indvI][1]) or last; 
			$eP_2   = $blkH->{'blk'}[$i][2]; 
			$blkJ_2e= $i; 
			$eP_2 - $eP_0 <= $glob{'len_dbX'} or last; # No need to count the following, because there will not be any double crossing. 
		} elsif ( $blkH->{'blk'}[$i][1] - $eP_0 <= $glob{'len_dbX'} ) {
			&_is_same_genotype('u', $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
			$abh_2  = $blkH->{'blk'}[$i][5][$indvI][1]; 
			$sP_2   = $blkH->{'blk'}[$i][1]; 
			$eP_2   = $blkH->{'blk'}[$i][2]; 
			$blkJ_2s= $i; 
			$blkJ_2e= $i; 
		} else {
			last; 
		}
	}
	unless (defined $abh_2 and $eP_2 - $eP_0 <= $glob{'len_dbX'}) {
		# There is no double crossing found from given $blkJ_1, so no need to change; 
		return(0, $runs_cnt); # ( 1-changed|0-not_changed, $runs_iterated)
	}
	@blkJ_idx_2 = ($blkJ_2s .. $blkJ_2e); 
	$n2 = &get_blkIndvN2( $blkH, $indvI, \@blkJ_idx_2, $abh_2 ); 
	
	# Extend to the right for n3 from n2; 
	my ($sP_3, $eP_3, $blkJ_3s, $blkJ_3e); 
	for (my $i=$blkJ_2e+1; $i<@{$blkH->{'blk'}}; $i++) {
		$chID_1 eq $blkH->{'blk'}[$i][0] or last; 
		if ( defined $abh_3 ) {
			$blkH->{'blk'}[$i][1] - $eP_3 <= $glob{'len_dbX'} or last; 
			&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
			&_is_same_genotype($abh_3, $blkH->{'blk'}[$i][5][$indvI][1]) or last; 
			$eP_3   = $blkH->{'blk'}[$i][2]; 
			$blkJ_3e= $i; 
		} elsif ( $blkH->{'blk'}[$i][1] - $eP_0 <= $glob{'len_dbX'} ) {
			&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
			$sP_3   = $blkH->{'blk'}[$i][1]; 
			$eP_3   = $blkH->{'blk'}[$i][2]; 
			$abh_3  = $blkH->{'blk'}[$i][5][$indvI][1]; 
			$blkJ_3s= $i; 
			$blkJ_3e= $i; 
		} else {
			last; 
		} 
	}
	unless (defined $abh_3) {
		# No n3 segment found. 
		return(0, $runs_cnt); 
	}
	@blkJ_idx_3 = ($blkJ_3s .. $blkJ_3e); 
	$n3 = &get_blkIndvN2( $blkH, $indvI, \@blkJ_idx_3, $abh_3 ); 
	
	my @has_changed; 
	
	# Check when ( $abh_1 eq $abh_3 ); 
#			if ( $blkH->{'offs'}[$indvI] eq 'F15-4231-rep1' and $blkH->{'blk'}[$blkJ_idx_1[0]][1] <= 36345) {
#				warn "S1: $blkH->{'blk'}[$blkJ_idx_1[0]][1] - $blkH->{'blk'}[$blkJ_idx_1[-1]][2]; @blkJ_idx_1\n"; 
#				warn "S2: $blkH->{'blk'}[$blkJ_idx_2[0]][1] - $blkH->{'blk'}[$blkJ_idx_2[-1]][2]; @blkJ_idx_2\n"; 
#				warn "S3: $blkH->{'blk'}[$blkJ_idx_3[0]][1] - $blkH->{'blk'}[$blkJ_idx_3[-1]][2]; @blkJ_idx_3\n"; 
#				warn "$n1 $n2 $n3\n$abh_1 $abh_2 $abh_3\n"; 
#			}
	if ( &_is_same_genotype($abh_1, $abh_3) ) {
		if      ( &_compare_siteN($n1, $n2) == 1 and &_compare_siteN($n3, $n2) == 1 ) {
			push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, $abh_1, $n2 ) ); 
		} elsif ( &_compare_siteN($n1, $n2) == 1 and &_compare_siteN($n3, $n2) == 0 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'u', $n2 ) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'h', $n2 ) ); 
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'h', $n3 ) ); 
			}
		} elsif ( &_compare_siteN($n3, $n2) == 1 and &_compare_siteN($n1, $n2) == 0 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'u', $n2 ) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'h', $n2 ) ); 
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'h', $n1 ) ); 
			}
		} elsif ( &_compare_siteN($n1, $n2) == 1 and &_compare_siteN($n3, $n2) == -1 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, $abh_1, $n2 ) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'u', $n3 ) ); 
			}
		} elsif ( &_compare_siteN($n3, $n2) == 1 and &_compare_siteN($n1, $n2) == -1 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, $abh_1, $n2 ) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'u', $n1 ) ); 
			}
		} elsif ( &_compare_siteN($n1, $n2) == 0 and &_compare_siteN($n3, $n2) == 0 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, $abh_1, $n2) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'h', $n1) ); 
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'h', $n2) ); 
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'h', $n3) ); 
			}
		} elsif ( &_compare_siteN($n1, $n2) == 0 and &_compare_siteN($n3, $n2) == -1 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, $abh_1, $n2) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'u', $n3) ); 
			}
		} elsif ( &_compare_siteN($n3, $n2) == 0 and &_compare_siteN($n1, $n2) == -1 ) {
			if (&_is_same_genotype($abh_2, 'u')) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, $abh_1, $n2) ); 
			} else {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'u', $n1) ); 
			}
		} elsif ( &_compare_siteN($n1, $n2) == -1 and &_compare_siteN($n3, $n2) == -1 ) {
			unless ( &_is_same_genotype($abh_2, 'u') ) {
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'u', $n1) ); 
				push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'u', $n3) ); 
			}
		} else {
			my $x1 = &_compare_siteN($n1, $n2); 
			my $x2 = &_compare_siteN($n3, $n2); 
			&stopErr("[Err] Unknown [n1($n1) cmp n2($n2) = $x1] and [n3($n3) cmp n2($n2) = $x2]\n"); 
		}
		return(&chk_change(@has_changed), $runs_cnt); 
	} else {
		# $abh_1 (S1) is different from $abh_3 (S3/SX); 
		# So go on to check if we can find n4; 
	}
	
	# Extend to the right for n4 from n3; 
	&_is_same_genotype($abh_2, 'u') and return(0, $runs_cnt); 
	&_is_same_genotype($abh_3, 'u') and return(0, $runs_cnt); 
	unless ($eP_3 - $eP_0 <= $glob{'len_dbX'}) {
		return(0, $runs_cnt); 
	}
	my ($sP_4, $eP_4, $blkJ_4s, $blkJ_4e); 
	for (my $i=$blkJ_3e+3; $i<@{$blkH->{'blk'}}; $i++) {
		$chID_1 eq $blkH->{'blk'}[$i][0] or last; 
		if ( defined $abh_4 ) {
			$blkH->{'blk'}[$i][1] - $eP_4 <= $glob{'len_dbX'} or last; 
			&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
			&_is_same_genotype($abh_4, $blkH->{'blk'}[$i][5][$indvI][1]) or last; 
			$eP_4   = $blkH->{'blk'}[$i][2]; 
			$blkJ_4e= $i; 
		} elsif ( $blkH->{'blk'}[$i][1] - $eP_0 <= $glob{'len_dbX'} ) {
			&_is_same_genotype('u',    $blkH->{'blk'}[$i][5][$indvI][1]) and next; 
			$sP_4   = $blkH->{'blk'}[$i][1]; 
			$eP_4   = $blkH->{'blk'}[$i][2]; 
			$abh_4  = $blkH->{'blk'}[$i][5][$indvI][1]; 
			$blkJ_4s= $i; 
			$blkJ_4e= $i; 
		} else {
			last; 
		}
	}
	unless (defined $abh_4) {
		return(0, $runs_cnt); 
	}
	@blkJ_idx_4 = ($blkJ_4s .. $blkJ_4e); 
	&_is_same_genotype($abh_2, $abh_4) and return(0, $runs_cnt); # Goto compare site_abh_2 to site_abh_4; 
	$n4 = &get_blkIndvN2( $blkH, $indvI, \@blkJ_idx_4, $abh_4 ); 
	
	# Check when ( $abh_2 ne $abh_4 )
	if      ( &_compare_siteN($n1,$n2) == 1 and &_compare_siteN($n2, $n3) == 0 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'h', $n2) ); 
	} elsif ( &_compare_siteN($n1,$n2) == 1 and &_compare_siteN($n2, $n3) == 1 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'u', $n3) ); 
	} elsif ( &_compare_siteN($n1,$n2) == 1 and &_compare_siteN($n2, $n3) == -1 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'u', $n2) ); 
	} elsif ( &_compare_siteN($n1,$n2) == 0 and &_compare_siteN($n2, $n3) == 1 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'u', $n3) ); 
	} elsif ( &_compare_siteN($n1,$n2) == 0 and &_compare_siteN($n2, $n3) == 0 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'h', $n1) ); 
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'h', $n2) ); 
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'h', $n3) ); 
	} elsif ( &_compare_siteN($n1,$n2) == 0 and &_compare_siteN($n2, $n3) == -1 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'h', $n1) ); 
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_2, 'h', $n2) ); 
	} elsif ( &_compare_siteN($n1,$n2) == -1 and &_compare_siteN($n2, $n3) == 0 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'u', $n1) ); 
	} elsif ( &_compare_siteN($n1,$n2) == -1 and &_compare_siteN($n2, $n3) == 1 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'u', $n1) ); 
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_3, 'u', $n3) ); 
	} elsif ( &_compare_siteN($n1,$n2) == -1 and &_compare_siteN($n2, $n3) == -1 ) {
		push( @has_changed, &replace_abh( $blkH, $indvI, \@blkJ_idx_1, 'u', $n1) ); 
	} else {
		return( &chk_change(@has_changed), $runs_cnt ); 
	}
	return( &chk_change(@has_changed), $runs_cnt ); 
}# smooth_blkJ() 

sub chk_change {
	my $is_changed = 0; 
	for my $t (@_) {
		$t == 0 and next; 
		$is_changed = $t; 
		last; 
	}
	return($is_changed); 
}# chk_change() 

sub replace_abh {
	my ($blkH, $indvI, $blkJ_ar, $new_abh, $blkJ_num) = @_; 
	$blkJ_num //= 0; 
	$blkJ_num >= $glob{'fixGeno_siteN'} and return(0); 
	my $has_changed = 0; 
	for my $i (@$blkJ_ar) {
		&_is_same_genotype( $blkH->{'blk'}[$i][5][$indvI][1], $new_abh ) and next; 
		$has_changed = 1; 
		$blkH->{'blk'}[$i][5][$indvI][1] = $new_abh; 
		$blkH->{'blk'}[$i][5][$indvI][2] = &_is_same_genotype($new_abh, $blkH->{'blk'}[$i][5][$indvI][0]) ? $blkH->{'blk'}[$i][4] : 0 ; 
	}
	return($has_changed); 
}# replace_abh() 


sub get_abh {
	&chk_para('get_abh'); 
	my $inTabH = &openFH( $glob{'vcf_tab'}, '<' ); 
	my $need_header = 1; 
	my $tab_header; 
	my @offs_colN; 
	if ( defined $glob{'offs_colN'} ) {
		@offs_colN = &mathSunhh::_parseCol( $glob{'offs_colN'} ); 
		@offs_colN > 0 or &stopErr("[Err] Input -offs_colN [$glob{'offs_colN'}] doesn't work.\n"); 
	}
	LINE: 
	while (readline($inTabH)) {
		chomp($_); 
		my @ta = &splitL("\t", $_); 
		if ($ta[0] =~ m!^(#CHROM|chr|chromosome|chrID|CHROM)$!i) {
			# This is a header of file. 
			$tab_header //= [ @ta ]; 
			next; 
		}
		unless ( @offs_colN > 0 ) {
			for (my $i=3; $i<@ta; $i++) {
				$i == $glob{'P1_colN'} and next; 
				$i == $glob{'P2_colN'} and next; 
				push(@offs_colN, $i); 
			}
		}
		# Get P1 and P2 allele; 
		my @p1_al = &SNP_tbl::tab_allele( $ta[$glob{'P1_colN'}] ); 
		my @p2_al = &SNP_tbl::tab_allele( $ta[$glob{'P2_colN'}] ); 
		if ( $glob{'infer_PG'} ) {
			@p1_al == 1 or do { &tsmsg("[Rec] Skip line with bad P1 genotype [$ta[$glob{'P1_colN'}]]: $_\n"); next; }; 
			@p2_al == 1 or do { &tsmsg("[Rec] Skip line with bad P2 genotype [$ta[$glob{'P2_colN'}]]: $_\n"); next; }; 
			( $p1_al[0][0] eq '.' and $p2_al[0][0] eq '.' ) and do { &tsmsg("[Rec] Skip line without parent allele information: $_\n"); next; }; 
		} else {
			(@p1_al == 1 and $p1_al[0][0] ne '.') or do { &tsmsg("[Rec] Skip line with bad P1 genotype [$ta[$glob{'P1_colN'}]]: $_\n"); next; }; 
			(@p2_al == 1 and $p2_al[0][0] ne '.') or do { &tsmsg("[Rec] Skip line with bad P2 genotype [$ta[$glob{'P2_colN'}]]: $_\n"); next; }; 
			$p1_al[0][0] ne $p2_al[0][0] or do { &tsmsg("[Rec] Skip line with same parent genotypes [$p1_al[0][0] VS. $p2_al[0][0]]: $_\n"); next; }; 
		}
		$p1_al[0][0] eq $p2_al[0][0] and do { &tsmsg("[Rec] Skip line with same parent alleles [$p1_al[0][0] VS. $p2_al[0][0]]: $_\n"); next LINE; }; # If both parents are './.', I'll skip this line too. 
		# Get offsprings' genotypes; 
		my %off_geno; 
		for my $i (@offs_colN) {
			$off_geno{'cntN'}{'offs_totalN'} ++; 
			$off_geno{'offs_al'}[$i] = [ &SNP_tbl::tab_allele( $ta[$i] ) ]; 
			for my $tr1 (@{$off_geno{'offs_al'}[$i]}) {
				$tr1->[0] eq '.' and next; 
				$off_geno{'al2cnt'}{ $tr1->[0] } += $tr1->[1]; 
			}
		}
		my @off_al_arr = sort { $off_geno{'al2cnt'}{$b} <=> $off_geno{'al2cnt'}{$a} || $a cmp $b } keys %{ $off_geno{'al2cnt'} }; 
		if ( @off_al_arr == 0 or @off_al_arr > 2 ) {
			&tsmsg("[Wrn] Skip line with abnormal alleles [$p1_al[0][0] VS. $p2_al[0][0]] [@off_al_arr] : $_\n"); 
			next LINE; 
		}
		# Infer parent alleles if required; 
		if ( $glob{'infer_PG'} ) {
			if ( $p1_al[0][0] eq '.' and $p2_al[0][0] ne '.' ) {
				for my $t_al ( @off_al_arr ) {
					$t_al eq $p2_al[0][0] and next; 
					$p1_al[0][0] = $t_al; 
					last; 
				}
			} elsif ( $p2_al[0][0] eq '.' and $p1_al[0][0] ne '.' ) {
				for my $t_al ( @off_al_arr ) {
					$t_al eq $p1_al[0][0] and next; 
					$p2_al[0][0] = $t_al; 
					last; 
				}
			}
		}
		( $p1_al[0][0] ne '.' and $p2_al[0][0] ne '.' ) or do { &tsmsg("[Rec] Skip line with missing parent alleles [$p1_al[0][0] VS. $p2_al[0][0]]: $_\n"); next LINE; }; 
		$off_geno{'pp_al_class'} = &SNP_tbl::tab_class_PP_al( \@p1_al, \@p2_al ); 
		# Check if offsprings' genotypes are consistent with parents. 
		for my $t_al ( @off_al_arr ) {
			if ( $t_al eq $p1_al[0][0] ) {
				$off_geno{'al2parent'}{$t_al} = 'P1'; 
			} elsif ( $t_al eq $p2_al[0][0] ) {
				$off_geno{'al2parent'}{$t_al} = 'P2'; 
			} elsif ( !$glob{'keep_abnormal'} ) { 
				&tsmsg("[Rec] Skip line with inconsistent parent/offsprings alleles [$p1_al[0][0] $p2_al[0][0] | @off_al_arr]: $_\n"); 
				next LINE;
			}
		}# for my $t_al 
		# Output tab_header; 
		if ( $need_header ) {
			$tab_header //= [ qw/chr pos ref/ ]; 
			$tab_header->[$glob{'P1_colN'}] //= 'P1'; $tab_header->[$glob{'P1_colN'}] eq '' and $tab_header->[$glob{'P1_colN'}] = 'P1'; 
			$tab_header->[$glob{'P2_colN'}] //= 'P2'; $tab_header->[$glob{'P2_colN'}] eq '' and $tab_header->[$glob{'P2_colN'}] = 'P2'; 
			for my $i ( @offs_colN ) {
				$tab_header->[$i] //= "offs_c$i"; $tab_header->[$i] eq '' and $tab_header->[$i] = "offs_c$i"; 
			}
			print STDOUT join("\t", @{$tab_header}[0,1,2,$glob{'P1_colN'}, $glob{'P2_colN'}, @offs_colN])."\n"; 
			$need_header = 0; 
		}
		# Output tab_sites; 
		for my $i ( @offs_colN ) {
			$off_geno{'offs_class'}[$i] = &SNP_tbl::tab_class_off_al( $off_geno{'pp_al_class'}, $off_geno{'offs_al'}[$i] ); 
			if (      $off_geno{'offs_class'}[$i] eq 'miss' ) {
				$off_geno{'offs_abh'}[$i] = '-'; 
			} elsif ( $off_geno{'offs_class'}[$i] eq 'homo_P1_parent' ) {
				$off_geno{'offs_abh'}[$i] = 'a'; 
			} elsif ( $off_geno{'offs_class'}[$i] eq 'homo_P2_parent' ) {
				$off_geno{'offs_abh'}[$i] = 'b'; 
			} elsif ( $off_geno{'offs_class'}[$i] eq 'hete_both_parent' ) {
				$off_geno{'offs_abh'}[$i] = 'h'; 
			} else {
				&stopErr("[Err] Bad allele [$off_geno{'offs_al'}[$i]] type [$off_geno{'offs_class'}[$i]] in line: $_\n"); 
			}
		}# for my $i ( @offs_colN )
		print STDOUT join("\t", 
			@ta[0,1,2], 
			&SNP_tbl::tab_allele_to_genotype(\@p1_al), 
			&SNP_tbl::tab_allele_to_genotype(\@p2_al), 
			@{$off_geno{'offs_abh'}}[@offs_colN]
		)."\n"; 
	}# End while (readline($inTabH)) 
	close($inTabH); 
	return; 
}# get_abh() 


### Format of vcf_tab file : 
##### #CHROM          POS     REF     F15-3866-rep1   F15-3866-rep2   F15-3866-rep3   F15-3882-rep1   F15-3882-rep2   F15-3882-rep3   F15-3886-rep1
##### SL2.40ch00      1147773 T       T/T             T/T             T/T             T/T             ./.     T/T     T/T             T/T
##### SL2.40ch00      1149665 T       T/T             T/T             T/T             T/T             T/T     T/T     T/T             T/T
##### SL2.40ch00      1149709 A       A/A             A/A             A/A             A/A             A/A     A/A     A/A             A/A
##### SL2.40ch01      1894598 TC      TC/TC           TC/TC           TC/TC           TC/TC           TC/TC   TC/TC   TC/TC           TC/TC
##### SL2.40ch01     75102180 CAA     CA/CA           C/CA            CAA/CA          CA/CA           CA/CA   CAA/CA  CA/CA           C/CA
sub cnt_pop_allele {
	&chk_para('cnt_pop_allele'); 
	my $inTabH = &openFH( $glob{'vcf_tab'}, '<' ); 
	my $need_header = 1; 
	my $tab_header; 
	my @offs_colN; 
	if ( defined $glob{'offs_colN'} ) {
		@offs_colN = &mathSunhh::_parseCol( $glob{'offs_colN'} ); 
		@offs_colN > 0 or &stopErr("[Err] Input -offs_colN [$glob{'offs_colN'}] doesn't work.\n"); 
	}
	LINE: 
	while (readline($inTabH)) {
		chomp($_); 
		my @ta = &splitL("\t", $_); 
		if ($ta[0] =~ m!^(#CHROM|chr|chromosome|chrID|CHROM)$!i) {
			# This is a header of file. 
			$tab_header //= [ @ta ]; 
			next; 
		}
		unless ( @offs_colN > 0 ) {
			for (my $i=3; $i<@ta; $i++) {
				$i == $glob{'P1_colN'} and next; 
				$i == $glob{'P2_colN'} and next; 
				push(@offs_colN, $i); 
			}
		}
		# Get P1 and P2 allele; 
		my @p1_al = &SNP_tbl::tab_allele( $ta[$glob{'P1_colN'}] ); 
		my @p2_al = &SNP_tbl::tab_allele( $ta[$glob{'P2_colN'}] ); 
		unless ( $glob{'keep_abnormal'} ) {
			if ( $glob{'infer_PG'} ) {
				@p1_al == 1 or do { &tsmsg("[Rec] Skip line with bad P1 genotype [$ta[$glob{'P1_colN'}]]: $_\n"); next; }; 
				@p2_al == 1 or do { &tsmsg("[Rec] Skip line with bad P2 genotype [$ta[$glob{'P2_colN'}]]: $_\n"); next; }; 
				( $p1_al[0][0] eq '.' and $p2_al[0][0] eq '.' ) and do { &tsmsg("[Rec] Skip line without parent allele information: $_\n"); next; }; 
			} else {
				(@p1_al == 1 and $p1_al[0][0] ne '.') or do { &tsmsg("[Rec] Skip line with bad P1 genotype [$ta[$glob{'P1_colN'}]]: $_\n"); next; }; 
				(@p2_al == 1 and $p2_al[0][0] ne '.') or do { &tsmsg("[Rec] Skip line with bad P2 genotype [$ta[$glob{'P2_colN'}]]: $_\n"); next; }; 
				$p1_al[0][0] ne $p2_al[0][0] or do { &tsmsg("[Rec] Skip line with same parent genotypes [$p1_al[0][0] VS. $p2_al[0][0]]: $_\n"); next; }; 
			}
		}
		# Get offsprings' genotypes; 
		my %off_geno; 
		for my $i (@offs_colN) {
			$off_geno{'cntN'}{'offs_totalN'} ++; 
			$off_geno{'offs_al'}[$i] = [ &SNP_tbl::tab_allele( $ta[$i] ) ]; 
			for my $tr1 (@{$off_geno{'offs_al'}[$i]}) {
				$tr1->[0] eq '.' and next; 
				$off_geno{'al2cnt'}{ $tr1->[0] } += $tr1->[1]; 
			}
		}
		my @off_al_arr = sort { $off_geno{'al2cnt'}{$b} <=> $off_geno{'al2cnt'}{$a} || $a cmp $b } keys %{ $off_geno{'al2cnt'} }; 
		if ( !$glob{'keep_abnormal'} and ( @off_al_arr == 0 or @off_al_arr > 2 ) ) {
			&tsmsg("[Wrn] Skip line with abnormal alleles [$p1_al[0][0] VS. $p2_al[0][0]] [@off_al_arr] : $_\n"); 
			next LINE; 
		}
		# Infer parent alleles if required; 
		if ( $glob{'infer_PG'} ) {
			if ( $p1_al[0][0] eq '.' and $p2_al[0][0] ne '.' ) {
				for my $t_al ( @off_al_arr ) {
					$t_al eq $p2_al[0][0] and next; 
					$p1_al[0][0] = $t_al; 
					last; 
				}
			} elsif ( $p2_al[0][0] eq '.' and $p1_al[0][0] ne '.' ) {
				for my $t_al ( @off_al_arr ) {
					$t_al eq $p1_al[0][0] and next; 
					$p2_al[0][0] = $t_al; 
					last; 
				}
			}
			unless ( $glob{'keep_abnormal'} ) {
				( $p1_al[0][0] ne '.' and $p2_al[0][0] ne '.' ) or do { &tsmsg("[Rec] Skip line with inconsistent parent/offsprings alleles [$p1_al[0][0] $p2_al[0][0] | @off_al_arr]: $_\n"); next LINE; }; 
			}
		}
		$off_geno{'pp_al_class'} = &SNP_tbl::tab_class_PP_al( \@p1_al, \@p2_al ); 
		# Check if offsprings' genotypes are consistent with parents. 
		for my $t_al ( @off_al_arr ) {
			if ( $t_al eq $p1_al[0][0] ) {
				$off_geno{'al2parent'}{$t_al} = 'P1'; 
			} elsif ( $t_al eq $p2_al[0][0] ) {
				$off_geno{'al2parent'}{$t_al} = 'P2'; 
			} elsif ( !$glob{'keep_abnormal'} ) { 
				&tsmsg("[Rec] Skip line with inconsistent parent/offsprings alleles [$p1_al[0][0] $p2_al[0][0] | @off_al_arr]: $_\n"); 
				next LINE;
			}
		}# for my $t_al 
		# Count offsprings alleles; 
		for my $i (@offs_colN) {
			if ( $#{$off_geno{'offs_al'}[$i]} == 0 ) {
				$off_geno{'geno_cntN'}{"$off_geno{'offs_al'}[$i][0][0]/$off_geno{'offs_al'}[$i][0][0]"} ++; 
			} else {
				$off_geno{'geno_cntN'}{"$off_geno{'offs_al'}[$i][0][0]/$off_geno{'offs_al'}[$i][1][0]"} ++; 
			}
			$off_geno{'offs_class'}[$i] = &SNP_tbl::tab_class_off_al( $off_geno{'pp_al_class'}, $off_geno{'offs_al'}[$i] ); 
			if (      $off_geno{'offs_class'}[$i] eq 'miss' ) {
				$off_geno{'cntN'}{'miss_N'} ++; 
			} elsif ( $off_geno{'offs_class'}[$i] eq 'homo_P1_parent' ) {
				$off_geno{'cntN'}{'homo_P1'} ++; 
			} elsif ( $off_geno{'offs_class'}[$i] eq 'homo_P2_parent' ) {
				$off_geno{'cntN'}{'homo_P2'} ++; 
			} elsif ( $off_geno{'offs_class'}[$i] eq 'hete_both_parent' ) {
				$off_geno{'cntN'}{'hete_PP'} ++; 
			} elsif ( $off_geno{'offs_class'}[$i] =~ m!_non_parent$! ) {
				$off_geno{'cntN'}{'non_PP'} ++; 
			} elsif ( $off_geno{'offs_class'}[$i] =~ m!_bad_parent$! ) {
				$off_geno{'cntN'}{'bad_PP'} ++; 
			} else {
				$off_geno{'cntN'}{'other'}{ $off_geno{'offs_class'}[$i] } ++; 
			}
		}
		$off_geno{'cntN'}{'other'}{'geno_cntInOffs'} = join(':', map { "|$_|=$off_geno{'geno_cntN'}{$_}" } sort { $off_geno{'geno_cntN'}{$b} <=> $off_geno{'geno_cntN'}{$a} || $a cmp $b } keys %{$off_geno{'geno_cntN'}}); 
		
		# Output header of file. 
		if ( $need_header ) {
			$tab_header //= [ 'chr', 'pos', 'ref' ]; 
			$tab_header->[$glob{'P1_colN'}] //= 'P1'; $tab_header->[$glob{'P1_colN'}] eq '' and $tab_header->[$glob{'P1_colN'}] = 'P1'; 
			$tab_header->[$glob{'P2_colN'}] //= 'P2'; $tab_header->[$glob{'P2_colN'}] eq '' and $tab_header->[$glob{'P2_colN'}] = 'P2'; 
			print STDOUT join("\t", @{$tab_header}[0,1,2, $glob{'P1_colN'}, $glob{'P2_colN'}], qw/offs_N miss homo_P1 homo_P2 hete_PP non_PP bad_PP others/)."\n"; 
			$need_header = 0; 
		}
		for my $tk (qw/offs_totalN miss_N homo_P1 homo_P2 hete_PP non_PP bad_PP/) {
			$off_geno{'cntN'}{$tk} //= 0; 
		}
		$off_geno{'other'} //= { 'NA'=>'NA' }; 
		print STDOUT join("\t", 
			@ta[0,1,2,$glob{'P1_colN'}, $glob{'P2_colN'}], 
			@{$off_geno{'cntN'}}{qw/offs_totalN miss_N homo_P1 homo_P2 hete_PP non_PP bad_PP/}, 
			join(";;", map { "$_=$off_geno{'cntN'}{'other'}{$_}" } sort keys %{$off_geno{'cntN'}{'other'}})
		) . "\n"; 
		
	}# End while 
	close($inTabH); 
	return; 
}# cnt_pop_allele() 


####################################################################################################
#       Inner sub-routines 
####################################################################################################
sub initial_glob {
	$glob{'P1_colN'} = 3; 
	$glob{'P2_colN'} = 4; 
	$glob{'task'}    = 'NA'; 
	{
		# For task:blk_abh; 
		$glob{'len_jnBlk'}     = 1000000; 
		$glob{'skipU_jnBlk'}   = undef(); 
		# For -smooth_geno in task:blk_abh; 
		$glob{'smooth_round'}  = 10; 
		$glob{'len_noX'}       = 1000; 
		$glob{'len_singleP'}   = 200; 
		$glob{'len_dbX'}       = 100000; 
		$glob{'fixGeno_siteN'} = 10; 
		$glob{'bigger_copyN'}  = 1.5; 
		$glob{'bigger_miniN'}  = 3; 
		$glob{'bigger_diffN'}  = 3; 
		$glob{'window_hRatio'} = 0.2; 
		$glob{'window_aRatio'} = 0.3; 
		$glob{'window_bRatio'} = 0.3; 
	}
	{
		# For task:rejoin_abhblk
		$glob{'min_isiteN2'}   = 10; 
	}
}# setupt_glob() 


sub chk_para {
	my ($task) = @_; 
	if      ( $task eq 'cnt_pop_allele') {
		for my $tk (qw/vcf_tab P1_colN P2_colN/) {
			defined $glob{$tk} or &usage("[Err] Need to assign [$tk\]n"); 
		}
	} elsif ( $task eq 'get_abh' ) { 
		for my $tk (qw/vcf_tab P1_colN P2_colN/) {
			defined $glob{$tk} or &usage("[Err] Need to assign [$tk]\n"); 
		}
		$glob{'keep_abnormal'} and &tsmsg("[Wrn] -keep_abnormal is not valid for -task get_abh\n"); 
	} elsif ( $task eq 'blk_abh' ) { 
		defined $glob{'abh_tab'} or &usage("[Err] Require -abh_tab file for -task blk_abh.\n"); 
		if ( $glob{'smooth_geno'} ) {
			$glob{'len_singleP'} > 0 or &usage("[Err] -len_singleP [$glob{'len_singleP'}] must be bigger than 0\n"); 
			for my $tk (qw/len_jnBlk len_noX len_singleP len_dbX fixGeno_siteN bigger_copyN bigger_miniN bigger_diffN/) {
				$glob{$tk} >= 0 or &usage("[Err] Bad input of -$tk [$glob{$tk}]\n"); 
			}
		}
	} elsif ( $task eq 'rejoin_abhblk' ) {
		defined $glob{'abhblk_tab'} or &tsmsg("[Err] Require -abhblk_tab file for -task rejoin_abhblk.\n"); 
		for my $tk (qw/min_isiteN2/) {
			$glob{$tk} >= 0 or &usage("[Err] Bad -$tk [$glob{$tk}]\n"); 
		}
		for my $tk (qw/len_jnBlk len_singleP/) {
			$glob{$tk} >  0 or &usage("[Err] Bad -$tk [$glob{$tk}]\n"); 
		}
	} else {
	}
}# chk_para() 

sub usage {
	&tsmsg(@_); 
my $help_txt = <<HH; 
####################################################################################################
# perl $0 
# -help 
# 
# -task     [string_of_task] 
#             Could be : 
#  task:cnt_pop_allele
#   -vcf_tab          [filename] Required. From vcf-to-tab < in.vcf > filename ; 
#   -P1_colN          [$glob{'P1_colN'}] Column number of parent P1; 
#   -P2_colN          [$glob{'P2_colN'}] Column number of parent P2; 
#   -offs_colN        [''] Column numbers of offsprings; Use all except parents if not given. 
#   -infer_PG         [Boolean] If this is given, try to infer one parent genotype from high frequence offspring alleles. 
#   -keep_abnormal    [Boolean] If this is given, all lines will be output ignoring if it is bad line. 
# 
#  task:get_abh 
#   Same parameters to task:cnt_pop_allele except 'keep_abnormal'; 
#   Return to STDOUT with format : chID       <tab> POS    <tab> REF_genotype <tab> P1_genotype_a <tab> P2_genotype_b     <tab> offspring_1_abh <tab> offspring_2_abh .... 
#                        example : #CHROM     <tab> POS    <tab> REF          <tab> VF36          <tab> Slycopersicoides2 <tab> F15-3866-rep1   <tab> F15-3866-rep2   .... 
#                        example : SL2.40ch00 <tab> 548583 <tab> AG           <tab> AG/AG         <tab> A/A               <tab> a               <tab> -               .... 
#     'a' - for P1 genotype; 
#     'b' - for P2 genotype; 
#     'h' - for heterozygous genotype; 
#     '-' - for unknown/missing genotype; 
# 
#  task:blk_abh 
#   -abh_tab          [pop.abh]
# 
#   Input .abh file from task:get_abh, and then output blocks for each individual; 
#   Output format : 
#   
#   -len_jnBlk        [$glob{'len_jnBlk'}] Any two same blocks within this length will be combined; 
#   -skipU_jnBlk      [Boolean] If given, genotype 'u' will be ignored when joining blocks. 
#   -window_geno      [Boolean] If given, I'll try to correct site-genotype according to window's genotype. 
#     -len_window     [$glob{'len_window'}] Length within which we test window genotype. 
#     -sN2_window     [$glob{'sN2_window'}] Individual siteN2 less than this number has window's genotype 'u'. 
#     -window_round   [$glob{'window_round'}] Maximum number of rounds to window genotypes.
#     -window_hRatio  [$glob{'window_hRatio'}] If having this ratio of heterozygous sites, this window is heterozygous; 
#     -window_aRatio  [$glob{'window_aRatio'}] If having this ratio of 'a' sites, this window has 'a' allele; 
#     -window_bRatio  [$glob{'window_bRatio'}] If having this ratio of 'b' sites, this window has 'b' allele; 
#                       If there are 'a' and 'b' alleles in window, this window is heterozygous; 
#                       Sometimes the reference is much similar to one of the parent, causing some mapping problem from 
#                       one parent to another. In this case, we can try to use unequal -window_aRatio and -window_bRatio
#                       values to give more possibility of the divergent parent origin alleles. 
#   -smooth_geno      [Boolean] If given, I'll try to check for double crossing and edit the genotypes of blocks. 
#     -smooth_round   [$glob{'smooth_round'}] Maximum number of rounds to smooth genotypes. 
#     -len_singleP    [$glob{'len_singleP'}] Length within which we count at most one SNP site (marker); 
#     -len_noX        [$glob{'len_noX'}] Length within which there should be no crossing; 
#     -len_dbX        [$glob{'len_dbX'}] Length within which we check double crossing problem; 
#     -fixGeno_siteN  [$glob{'fixGeno_siteN'}] If there are so many sites supporting this genotype, we don't change it at all (Except genotype 'u'); 
# 
#     -bigger_copyN   [$glob{'bigger_copyN'}] 
#     -bigger_miniN   [$glob{'bigger_miniN'}] 
#     -bigger_diffN   [$glob{'bigger_diffN'}]
#       Used to compare two siteN. siteN_1 is bigger than siteN_2 if all the followings three are true : 
#         siteN_1 >= bigger_copyN * siteN_2 
#         siteN_1 >= \$bigger_miniN 
#         siteN_1 >= \$bigger_diffN + siteN_2 
#       If any of the above is not satisfied, siteN_1 is similar to siteN_2; 
#
#  task:rejoin_abhblk
#   -abhblk_tab       [filename.abh.jn] Required. 
#   -min_isiteN2      [$glob{'min_isiteN2'}] Any block with less than this number of individual siteN2 are ignored. 
####################################################################################################

HH
	&LogInforSunhh::usage($help_txt); 
}# usage 

=head1 _compare_siteN( siteN_1, siteN_2 ) 

Return : ( $string )

  $string is one of [ 1, -1, 0 ]; 
	'1'  means siteN_1 > siteN_2 ; 
	'-1' means siteN_1 < siteN_2 ; 
	'0'  means siteN_1 ~ siteN_2 ; 

=cut
sub _compare_siteN {
	# ($_[0], $_[1]); 
	if ( $_[0] >= $_[1] ) {
		if    ( $_[0] >= $glob{'bigger_copyN'} * $_[1] and $_[0] >= $glob{'bigger_miniN'} and $_[0] >= $glob{'bigger_diffN'} + $_[1] ) {
			return(1); 
		} else {
			return(0); 
		}
	} elsif ( $_[1] >= $glob{'bigger_copyN'} * $_[0] and $_[1] >= $glob{'bigger_miniN'} and $_[1] >= $glob{'bigger_diffN'} + $_[0] ) {
		return(-1); 
	} else {
		return(0); 
	}
	&stopErr("[Err] Why here @_\n"); 
}# _compare_siteN() 

sub _is_same_genotype {
	lc($_[0]) eq lc($_[1]) and return 1; 
	return 0; 
}# _is_same_genotype() 

sub _is_same_neighbor_site {
	# ([indv_1_geno, indv_2_geno, ...], [indv_1_geno, indv_2_geno, ...], $geno_i)
	# Checking [0] which means $raw_genotype. 
	# a == a , a != u; 
	$_[2] //= 0; 
	my $is_same = 1; 
	for (my $i=0; $i<@{$_[0]}; $i++) {
		&_is_same_genotype( $_[0][$i][ $_[2] ] , $_[1][$i][ $_[2] ] ) or do { $is_same = 0; last; }; 
	}
	return($is_same); 
}# _is_same_neighbor_site() 


sub _load_abhblk {
	my ($fn) = @_; 
	
	my %back; 
	my $order = 0; 
	my $fh = &openFH($fn, '<'); 
	<$fh>; 
	while (<$fh>) {
		chomp; 
		my @ta = &splitL("\t", $_); 
		$back{$ta[0]}{'offs'}  //= [ $ta[0] ]; 
		defined $back{$ta[0]}{'order'} or do { $order ++; $back{$ta[0]}{'order'} = $order; }; 
		my @fake_sites; 
		if ( $ta[5] == 0 ) {
		} elsif ( $ta[5] == 1 ) {
			@fake_sites = int(($ta[2]+$ta[3])/2); 
		} elsif ( $ta[5] == 2 ) {
			@fake_sites = @ta[2,3]; 
		} else {
			my $delt = int( ($ta[3]-$ta[2]-1)/($ta[5]-1) ); 
			for (my $k=0; $k<$ta[5]-1; $k++) {
				push(@fake_sites, $ta[2]+$k*$delt); 
			}
			push(@fake_sites, $ta[3]); 
		}
		
		push(@{$back{$ta[0]}{'blk'}}, 
			[@ta[1..3,6,7], 
			[ [$ta[4], $ta[4], $ta[5]] ], 
			[ @fake_sites ]
			]
		); 
	}
	close($fh); 
	
	return(\%back); 
}# _load_abhblk() 

=head1 _load_abh( $abh_filename ) 

Required : $glob{qw/len_singleP len_noX/}

  Input abh file must be sorted by position and must have one line header information. 

Return : (\%blk_abh)

	$blk_abh{'blk'} = 
	[ 
		$chID, 
		$blk_startP, 
		$blk_endP, 
		$siteN1, 
		$siteN2, 
		[$geno_indv_1, $geno_indv_2, ...], 
		[$site_1_P,       $site_2_P,       ...], 
	]; 
	Where 
	$geno_indv_1 = [ $raw_geno_indv_1, $inferred_geno_indv_1, $inferred_geno_indv_1_siteN2 ]
	
	
	$blk_abh{'offs'} = [ $indv_1_ID, $indv_2_ID, ... ]; 
	$blk_ahb{'P1'}   = $P1_ID; 
	$blk_ahb{'P2'}   = $P2_ID; 
	$blk_ahb{'len_noX'}     = $glob{'len_noX'}; 
	$blk_ahb{'len_singleP'} = $glob{'len_singleP'}; 
	$blk_ahb{'site2geno'}   = { $chrID => { $position => [ $ref_geno, $P1_geno, $P2_geno ], ... }, ... }; 

=cut
sub _load_abh {
	my ($fn) = @_; 
	
	&tsmsg("[Msg] Start loading abh file $fn\n"); 
	
	my %back; 
	$back{'len_noX'}     = $glob{'len_noX'}; 
	$back{'len_singleP'} = $glob{'len_singleP'}; 
	my $fh = &openFH($fn, '<'); 
	my $hl = <$fh>; 
	chomp($hl); 
	my @ha = &splitL("\t", $hl); 
	$back{'P1'} = $ha[3]; 
	$back{'P2'} = $ha[4]; 
	$back{'offs'} = [@ha[5 .. $#ha]]; 
	while (<$fh>) {

$. % 10000 == 1 and &tsmsg("[Msg] $. lines.\n"); 
	
		chomp($_); 
		my @ta = &splitL("\t", $_); 
		my ($chID, $posi, $refG, $P1G, $P2G) = @ta[0,1,2,3,4]; 
		defined $back{'site2geno'}{$chID}{$posi} and do { &tsmsg("[Wrn] Skip repeated position: [$chID $posi].\n"); next; }; 
		$back{'site2geno'}{$chID}{$posi} = [ $refG, $P1G, $P2G ]; 
		
		### Check and formatted offsprings genotypes 
		for my $tb ( @ta[5 .. $#ta] ) {
			$tb =~ m!^[abhu\-.]$!i or &stopErr("[Err] Unknown offspring genotype [$tb]\n"); 
			( $tb eq '-' or $tb eq '.' ) and $tb = 'u'; 
			$tb = lc($tb); 
		}
		
		### Setup blocks by 'len_noX'
		if ( 
			defined $back{'blk'} 
			and $back{'blk'}[-1][0] eq $chID 
			and ($posi-$back{'blk'}[-1][2]) <= $back{'len_noX'}
			and &_is_same_neighbor_site( $back{'blk'}[-1][5], [ map { [ $_, $_, 0 ] } @ta[5..$#ta] ], 0 )
		) {
			$back{'blk'}[-1][2] = $posi; 
			$back{'blk'}[-1][3] ++; 
			push(@{$back{'blk'}[-1][6]}, $posi); 
		} else {
			my @tc = map { [ $_, $_, 0 ] } @ta[5 .. $#ta]; 
			push(@{$back{'blk'}}, 
				[ $chID, $posi, $posi, 1, 0, 
					[@tc], 
					[$posi]
				]
			); 
		}
	}
	close($fh); 
	### Update $siteN2 in blocks; 
	for my $tb (@{$back{'blk'}}) {
		my $sP = $tb->[1]; 
		my $cnt = 0; 
		for my $ts (@{$tb->[6]}) {
			$ts >= $sP or next; 
			$cnt ++; 
			while ( $ts >= $sP ) {
				$sP += $back{'len_singleP'}; 
			}
		}
		$tb->[4] = $cnt; 
		for my $td ( @{$tb->[5]} ) {
			$td->[2] = $cnt; 
		}
	}
	
	&tsmsg("[Msg] Finish loading abh file $fn\n"); 
	
	return(\%back); 
}# _load_abh() 

sub output_blkH_i {
	my ($blkH) = @_; 
	my @offs_IDs = @{$blkH->{'offs'}}; 
	@offs_IDs == 1 or &stopErr("[Err] There should be one and only one offspring ID. offs_ID=[@offs_IDs]\n"); 
	my $offs_ID = $offs_IDs[0]; 
	for my $tb (@{$blkH->{'blk'}}) {
		my ($chID, $sP, $eP, $siteN1, $siteN2) = @{$tb}[0,1,2,3,4]; 
		my @offs_geno = @{$tb->[5][0]}; 
		# Output format : offspring_ID <tab> chID <tab> startPos <tab> endPos <tab> offspring_abh <tab> offspring_siteN2 <tab> pop_siteN1 <tab> pop_siteN2
		unless (defined $glob{'has_header_output_blkH_i'} and $glob{'has_header_output_blkH_i'} == 1) {
			print STDOUT join("\t", qw/Offs_ID ChrID StartP EndP Offs_abh Offs_siteN2 Pop_siteN1 Pop_siteN2/)."\n"; 
			$glob{'has_header_output_blkH_i'} = 1; 
		}
		print STDOUT join("\t", $offs_ID, $chID, $sP, $eP, $offs_geno[1], $offs_geno[2], $siteN1, $siteN2)."\n"; 
	}
	return; 
}# output_blkH_i() 

sub _join_blk_i {
	my ($blkH, $indvI, $len, $skipU) = @_; 
	$skipU //= 0; 
	my %back; 
	$back{'offs'} = [$blkH->{'offs'}[$indvI]]; 
	for my $tk (qw/P1 P2 len_noX len_singleP site2geno/) {
		$back{$tk} = $blkH->{$tk}; 
	}
	my ($si, $ei); 
	my $ci = -1; 
	
	for my $tb (@{$blkH->{'blk'}}) {
		$ci ++; 
		my ($chID, $sP, $eP) = @{$tb}[0,1,2]; 
		my @offs_geno = @{ $tb->[5][$indvI] }; 
		$skipU and &_is_same_genotype( $offs_geno[1], 'u' ) and next; 
		if ( 
			defined $ei 
			and $blkH->{'blk'}[$ei][0] eq $chID 
			and ($sP-$blkH->{'blk'}[$ei][2]) <= $len 
			and &_is_same_genotype( $blkH->{'blk'}[$ei][5][$indvI][1], $offs_geno[1] )
		) {
			# extend blk_jn 
			$ei = $ci; 
		} else {
			if ( defined $ei ) {
				# Update previous block; 
				my $siteN1       = &mathSunhh::_sum( map { $_->[3] } @{$blkH->{'blk'}}[$si .. $ei] ); 
				my $siteN2       = &get_blkSiteN2( $blkH, [$si .. $ei] ); 
				my $siteN2_indv  = &get_blkIndvN2( $blkH, $indvI, [$si .. $ei], $blkH->{'blk'}[$si][5][$indvI][1], 0 ); 
				push(@{$back{'blk'}}, 
					[
						$blkH->{'blk'}[$si][0], 
						$blkH->{'blk'}[$si][1], 
						$blkH->{'blk'}[$ei][2], 
						$siteN1, 
						$siteN2, 
						[ [ $blkH->{'blk'}[$si][5][$indvI][0], $blkH->{'blk'}[$si][5][$indvI][1], $siteN2_indv ] ], 
						[ map { @{$_->[6]} } @{$blkH->{'blk'}}[$si .. $ei] ]
					]
				); 
				# $si = $ei = undef(); 
			}# if ( defined $ei ) 
			$si = $ei = $ci; 
		}
	}#End for my $tb (@{$blkH->{'blk'}})
	if ( defined $ei ) {
		# Update previous block; 
		my $siteN1       = &mathSunhh::_sum( map { $_->[3] } @{$blkH->{'blk'}}[$si .. $ei] ); 
		my $siteN2       = &get_blkSiteN2( $blkH, [$si .. $ei] ); 
		my $siteN2_indv  = &get_blkIndvN2( $blkH, $indvI, [$si .. $ei], $blkH->{'blk'}[$si][5][$indvI][1], 0 ); 
		push(@{$back{'blk'}}, 
			[
				$blkH->{'blk'}[$si][0], 
				$blkH->{'blk'}[$si][1], 
				$blkH->{'blk'}[$ei][2], 
				$siteN1, 
				$siteN2, 
				[ [ $blkH->{'blk'}[$si][5][$indvI][0], $blkH->{'blk'}[$si][5][$indvI][1], $siteN2_indv ] ], 
				[ map { @{$_->[6]} } @{$blkH->{'blk'}}[$si .. $ei] ]
			]
		); 
		$si = $ei = undef(); 
	}# if ( defined $ei ) 
	
	return(\%back); 
}# _join_blk_i() 

sub get_blkIndvN2 {
	my ($blkH, $indvI, $blk_idx_ar, $new_abh, $chk_genoI) = @_; 
	$chk_genoI //= 1; 
	my $n2 = 0; 
	my $sP; 
	for my $i (@$blk_idx_ar) {
		my $tb = $blkH->{'blk'}[$i]; 
		$sP //= $tb->[1]; 
		&_is_same_genotype( $blkH->{'blk'}[$i][5][$indvI][$chk_genoI], $new_abh ) or next; 
		for my $ts (@{$tb->[6]}) {
			$ts >= $sP and $n2 ++; 
			while ( $ts >= $sP ) {
				$sP += $glob{'len_singleP'}; 
			}
		}
	}
	return($n2); 
}# get_blkIndvN2() 

sub get_blkSiteN2 {
	my ($blkH, $blk_idx_ar) = @_; 
	my $n2 = 0; 
	my $sP; 
	for my $i (@$blk_idx_ar) {
		my $tb = $blkH->{'blk'}[$i]; 
		$sP //= $tb->[1]; 
		for my $ts (@{$tb->[6]}) {
			$ts >= $sP and $n2 ++; 
			while ( $ts >= $sP ) {
				$sP += $glob{'len_singleP'}; 
			}
		}
	}
	return($n2); 
}# get_blkSiteN2() 

sub normal_abh {
	return(lc($_[0])); 
}# 
sub ratio_abh {
	my ($wind_abh) = @_; 
	my %cnt; 
	for (@$wind_abh) {
		$cnt{ $_->[0] } ++; 
	}
	for (qw/a b h u/) {
		$cnt{$_} //= 0; 
	}
	$cnt{'total'} = $cnt{'a'}+$cnt{'b'}+$cnt{'h'}; 
	if ($cnt{'total'} <= 0) {
		$cnt{'a_r'} = $cnt{'b_r'} = $cnt{'h_r'} = -1; 
	} else {
		$cnt{'a_r'} = $cnt{'a'} / $cnt{'total'}; 
		$cnt{'b_r'} = $cnt{'b'} / $cnt{'total'}; 
		$cnt{'h_r'} = $cnt{'h'} / $cnt{'total'}; 
	}
	return(\%cnt); 
}# ratio_abh() 
sub type_abhR {
	my ($wind_abhR) = @_; 
	my %has; 
	$wind_abhR->{'total'} < $glob{'sN2_window'} and return('u'); # 
	$wind_abhR->{'h_r'} >= $glob{'window_hRatio'} and return('h'); 
	$wind_abhR->{'a_r'} + $wind_abhR->{'h_r'}/2 >= $glob{'window_aRatio'} and $has{'a'} = 1; 
	$wind_abhR->{'b_r'} + $wind_abhR->{'h_r'}/2 >= $glob{'window_bRatio'} and $has{'b'} = 1; 
	defined $has{'a'} and defined $has{'b'} and return('h'); 
	defined $has{'a'} and return('a'); 
	defined $has{'b'} and return('b'); 
	return('u'); 
}# type_abhR() 
sub refine_wind {
	my ($w) = @_; 
	my @wind_abh = @{$w->{'wind_abh'}}; 
	if ($w->{'lr'} eq 'left') {
		my $abh_prev; 
		my $stop_i; 
		if ( &_is_same_genotype($w->{'wind_type'}, 'h') ) {
			for (my $i=0; $i<@wind_abh; $i++) {
				unless ( &_is_same_genotype($wind_abh[$i][0], 'u') ) {
					$abh_prev //= $wind_abh[$i][0]; 
					&_is_same_genotype( $abh_prev, $wind_abh[$i][0] ) or last; 
				}
				$stop_i = $i; 
			}
			if ( defined $stop_i ) {
				@wind_abh = @wind_abh[($stop_i+1) .. $#wind_abh]; 
			}
			my $t_r = &ratio_abh( \@wind_abh ); 
			my $t_t = &type_abhR( $t_r ); 
			if ( &_is_same_genotype($t_t, 'a') or &_is_same_genotype($t_t, 'b') ) {
				$w->{'wind_abh'}  = \@wind_abh; 
				$w->{'wind_abhR'} = $t_r; 
				$w->{'wind_type'} = $t_t; 
			} elsif ( &_is_same_genotype($t_t, 'h') or &_is_same_genotype($t_t, 'u') ) {
				@wind_abh = @{$w->{'wind_abh'}}; 
				$abh_prev = $stop_i = undef(); 
				
				for (my $i=0; $i<@wind_abh; $i++) {
					&_is_same_genotype( $wind_abh[$i][0], $w->{'wind_type'} ) and last; 
					unless ( &_is_same_genotype($wind_abh[$i][0], 'u') ) {
						$abh_prev //= $wind_abh[$i][0]; 
						&_is_same_genotype( $abh_prev, $wind_abh[$i][0] ) or last; 
					}
					$stop_i = $i; 
				}
				if ( defined $stop_i ) {
					@wind_abh = @wind_abh[($stop_i+1) .. $#wind_abh]; 
				}
				my $t_r = &ratio_abh( \@wind_abh ); 
				my $t_t = &type_abhR( $t_r ); 
				$w->{'wind_abh'}  = \@wind_abh; 
				$w->{'wind_abhR'} = $t_r; 
				$w->{'wind_type'} = $t_t; 
			} else {
				&stopErr("[Err] Error here [$t_t].\n"); 
			}
		} else {
			for (my $i=0; $i<@wind_abh; $i++) {
				&_is_same_genotype( $w->{'wind_type'}, $wind_abh[$i][0] ) and last; 
				$stop_i = $i; 
			}
			if ( defined $stop_i ) {
				@wind_abh = @wind_abh[($stop_i+1) .. $#wind_abh]; 
			}
			my $t_r = &ratio_abh( \@wind_abh ); 
			my $t_t = &type_abhR( $t_r ); 
			$w->{'wind_abh'}   = \@wind_abh; 
			$w->{'wind_abhR'}  = $t_r; 
			$w->{'wind_type'}  = $t_t; 
		}
		if (@wind_abh == 0) {
			$w->{'wind_S_i'} = $w->{'wind_E_i'} = 'NA'; 
		} else {
			for my $ti ($w->{'wind_S_i'} .. $w->{'wind_E_i'}) {
				$w->{'blkH'}{'blk'}[$ti][2] < $w->{'wind_abh'}[0][1] and next; 
				$w->{'wind_S_i'} = $ti; 
				last; 
			}
		}
	} elsif ($w->{'lr'} eq 'right') {
		my $abh_prev; 
		my $stop_i; 
		if ( &_is_same_genotype($w->{'wind_type'}, 'h') ) {
			for (my $i=$#wind_abh; $i>=0; $i--) {
				unless ( &_is_same_genotype($wind_abh[$i][0],'u') ) {
					$abh_prev //= $wind_abh[$i][0]; 
					&_is_same_genotype( $abh_prev, $wind_abh[$i][0] ) or last; 
				}
				$stop_i = $i; 
			}
			if ( defined $stop_i ) {
				@wind_abh = @wind_abh[($stop_i+1) .. $#wind_abh]; 
			}
			my $t_r = &ratio_abh( \@wind_abh ); 
			my $t_t = &type_abhR( $t_r ); 
			if ( &_is_same_genotype($t_t, 'a') or &_is_same_genotype($t_t, 'b') ) {
				$w->{'wind_abh'}  = \@wind_abh; 
				$w->{'wind_abhR'} = $t_r; 
				$w->{'wind_type'} = $t_t; 
			} elsif ( &_is_same_genotype($t_t, 'h') or &_is_same_genotype($t_t, 'u') ) {
				@wind_abh = @{$w->{'wind_abh'}}; 
				$abh_prev = $stop_i = undef(); 
				for (my $i=$#wind_abh; $i>=0; $i--) {
					&_is_same_genotype( $wind_abh[$i][0], $w->{'wind_type'} ) and last; 
					unless ( &_is_same_genotype($wind_abh[$i][0],'u') ) {
						$abh_prev //= $wind_abh[$i][0]; 
						&_is_same_genotype( $abh_prev, $wind_abh[$i][0] ) or last; 
					}
					$stop_i = $i; 
				}
				if ( defined $stop_i ) {
					@wind_abh = @wind_abh[($stop_i+1) .. $#wind_abh]; 
				}
				my $t_r = &ratio_abh( \@wind_abh ); 
				my $t_t = &type_abhR( $t_r ); 
				$w->{'wind_abh'}  = \@wind_abh; 
				$w->{'wind_abhR'} = $t_r; 
				$w->{'wind_type'} = $t_t; 
			} else {
				&stopErr("[Err] Error here 2 [$t_t].\n"); 
			}
		} else {
			for (my $i=$#wind_abh; $i>=0; $i--) {
				&_is_same_genotype( $w->{'wind_type'}, $wind_abh[$i][0] ) and last; 
				$stop_i = $i; 
			}
			if ( defined $stop_i ) {
				@wind_abh = @wind_abh[($stop_i+1) .. $#wind_abh]; 
			}
			my $t_r = &ratio_abh( \@wind_abh ); 
			my $t_t = &type_abhR( $t_r ); 
			$w->{'wind_abh'}   = \@wind_abh; 
			$w->{'wind_abhR'}  = $t_r; 
			$w->{'wind_type'}  = $t_t; 
		}
		if (@wind_abh == 0) {
			$w->{'wind_S_i'} = $w->{'wind_E_i'} = 'NA'; 
		} else {
			for my $ti (reverse( $w->{'wind_S_i'} .. $w->{'wind_E_i'} )) {
				$w->{'blkH'}{'blk'}[$ti][1] > $w->{'wind_abh'}[-1][1] and next; 
				$w->{'wind_E_i'} = $ti; 
				last; 
			}
		}
	} else {
		&stopErr("[Err] die here [$w->{'lr'}].\n"); 
	}
	
	return; 
}# refine_wind() 

sub type_wind {
	my ($blkH, $indvI, $blkJ_0, $lr) = @_; 
	my $chID_0 = $blkH->{'blk'}[$blkJ_0][0]; 
	my %back; 
	$back{'lr'}   = $lr; 
	$back{'blkH'} = $blkH; 
	my $cP; 
	if ( $lr eq 'left' ) {
		my $posE_0 = $blkH->{'blk'}[$blkJ_0][2]; 
		$back{'wind_E_i'} = $blkJ_0; 
		$back{'wind_S_i'} = $blkJ_0; 
		$cP = $posE_0; 
		for (; $back{'wind_S_i'}>=0; $back{'wind_S_i'}--) {
			$blkH->{'blk'}[$back{'wind_S_i'}][0] eq $chID_0 or last; 
			$posE_0 - $blkH->{'blk'}[$back{'wind_S_i'}][2] <= $glob{'len_window'} or last; 
			my $abh_c = &normal_abh( $blkH->{'blk'}[$back{'wind_S_i'}][5][$indvI][1] ); 
			for my $tp (reverse( @{$blkH->{'blk'}[$back{'wind_S_i'}][6]} )) {
				$posE_0 - $tp <= $glob{'len_window'} or last; 
				$cP >= $tp or next; 
				push(@{$back{'wind_abh'}}, [$abh_c, $tp]); 
				while ($cP >= $tp) {
					$cP -= $glob{'len_singleP'}; 
				}
			}
		}# End for (; 
		@{$back{'wind_abh'}} = reverse(@{$back{'wind_abh'}}); 
		$back{'wind_abhR'}   = &ratio_abh( $back{'wind_abh'} ); 
		$back{'wind_type'}   = &type_abhR( $back{'wind_abhR'} ); 
		&refine_wind(\%back); 
	} elsif ( $lr eq 'right' ) {
		my $blkJ_1 = $blkJ_0; 
		my $posS_1 = $blkH->{'blk'}[$blkJ_1][1]; 
		$back{'wind_E_i'} = $blkJ_1; 
		$back{'wind_S_i'} = $blkJ_1; 
		$cP = $posS_1; 
		for (; $back{'wind_E_i'}<@{$blkH->{'blk'}}; $back{'wind_E_i'}++) {
			$blkH->{'blk'}[$back{'wind_E_i'}][0] eq $chID_0 or last; 
			$blkH->{'blk'}[$back{'wind_E_i'}][1] - $posS_1 <= $glob{'len_window'} or last; 
			my $abh_c = &normal_abh( $blkH->{'blk'}[$back{'wind_E_i'}][5][$indvI][1] ); 
			for my $tp ( @{$blkH->{'blk'}[$back{'wind_E_i'}][6]} ) {
				$tp - $posS_1 <= $glob{'len_window'} or last; 
				$cP <= $tp or next; 
				push(@{$back{'wind_abh'}}, [$abh_c, $tp]); 
				while ( $cP <= $tp ) {
					$cP += $glob{'len_singleP'}; 
				}
			}
		}# End for (; 
		$back{'wind_E_i'} > $#{$blkH->{'blk'}} and $back{'wind_E_i'}--; 
		$back{'wind_abhR'} = &ratio_abh( $back{'wind_abh'} ); 
		$back{'wind_type'} = &type_abhR( $back{'wind_abhR'} ); 
		&refine_wind(\%back); 
	} else {
		&stopErr("[Err] bad lr [$lr]\n"); 
	}
	return(\%back); 
}# type_wind() 

