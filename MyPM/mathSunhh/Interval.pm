package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 



#=head2 ovl_len( $start1, $end1, $start2, $end2 )
#
#Required: 
# ( $start1, $end1, $start2, $end2 )
#
#Function : 
# return overlapped length of to regions [$start1,$end1] and [$start2,$end2]
#
#=cut
#sub ovl_len {
#	my $self = shift; 
#	my ($s1, $e1, $s2, $e2) = @_; 
#	($s1, $e1) = sort {$a <=> $b} ($s1, $e1); 
#	($s2, $e2) = sort {$a <=> $b} ($s2, $e2); 
#	if ($e1 < $s2 or $s1 > $e2) {
#		return 0; 
#	} else {
#		return &mathSunhh::min($e1, $e2) - &mathSunhh::max($s1, $s2) + 1; 
#	}
#}
=head2 ovl_region( $start1, $end1, $start2, $end2 )

Required: ($start1, $end1, $start2, $end2)

Function: Check overlapping region

Return  : ( $overlapped_bp_length, [ ovl_start, ovl_end ] )

=cut
sub ovl_region {
	my ($s1, $e1, $s2, $e2) = @_; 
	$s1 > $e1 and ($s1, $e1) = ($e1, $s1); 
	$s2 > $e2 and ($s2, $e2) = ($e2, $s2); 
	if ( $e1 < $s2 or $s1 > $e2 ) {
		return (0, []); 
	} else {
		my $ee = &mathSunhh::min($e1, $e2); 
		my $ss = &mathSunhh::max($s1, $s2); 
		return ( $ee - $ss + 1, [$ss, $ee] ); 
	}
}# ovl_region(); 

=head2 compare_number_list( $sePair_list1, $sePair_list2, 
 'compare'=>'same/ovl/nonovl', 
 'sort'=>'0/single/pair'
)

Required : 

 $sePair_list1 : [ [s11,e11], [s12,e12], ... ]
 $sePair_list2 : [ [s21,e21], [s22,e22], ... ]

Return   : 

 Case 'compare'='same'   : 1/0 : 1 for yes_same, 0 for different. 
 Case 'compare'='ovl'    : ( $ovlLen_nonDup, $ovlCnt_mayDup, \@ovlLoc ) 
                         : @ovlLoc=([overlap_S1, overlap_E1], [overlap_S2, overlap_E2], 
 Case 'compare'='nonovl' : ( [@spec1], [@spec2] ) 
                         : @spec1=( [nonOvl_lis1_S1, nonOvl_lis1_E1], [nonOvl_lis1_S1, nonOvl_lis1_E1], ... )
                         : @spec2=( [nonOvl_lis2_S1, nonOvl_lis2_E1], [nonOvl_lis2_S1, nonOvl_lis2_E1], ... )

Function : 

 Case 'sort'='0'         : Use the input order in list1/list2; 
 Case 'sort'='single'    : Sort numbers from small to large one by one; 
 Case 'sort'='pair'      : Sort numbers from small to large as pairs, but no sorting within any pair. 

=cut
sub compare_number_list {
	my $list1 = shift; 
	my $list2 = shift; 
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	$parm{'compare'} = $parm{'compare'} // 'same'; # same/ovl/nonovl
	$parm{'compare'} = lc( $parm{'compare'} ); 
	$parm{'sort'} = $parm{'sort'} // 0; # 0/single/pair
	$parm{'sort'} = lc( $parm{'sort'} ); 
	my (@srt_lis1, @srt_lis2); 
	for ( @$list1 ) {
		if ( ref($_) eq "" ) {
			push(@srt_lis1, $_); 
		} elsif ( ref($_) eq 'ARRAY' ) {
			my @ta = @$_; 
			@ta > 2 and @ta = @ta[0,1]; 
			push(@srt_lis1, @ta); 
		} else {
		}
	}
	for ( @$list2 ) {
		if ( ref($_) eq "" ) {
			push(@srt_lis2, $_); 
		} elsif ( ref($_) eq 'ARRAY' ) {
			my @ta = @$_; 
			@ta > 2 and @ta = @ta[0,1]; 
			push(@srt_lis2, @ta); 
		} else {
		}
	}
	
	if ($parm{'sort'} =~ m/^s(ingle)?$/) {
		@srt_lis1 = map { [$_] } sort { $a <=> $b } @srt_lis1; 
		@srt_lis2 = map { [$_] } sort { $a <=> $b } @srt_lis2; 
	} elsif ( $parm{'sort'} =~ m/^p(air(ed)?)?$/) {
		scalar( @srt_lis1 ) % 2 == 0 or &stopErr("[Err] Input list1 is not even.\n"); 
		scalar( @srt_lis2 ) % 2 == 0 or &stopErr("[Err] Input list2 is not even.\n"); 
		my (@tmp_1, @tmp_2); 
		for (my $i=0; $i<@srt_lis1; $i+=2) { push(@tmp_1, [ @srt_lis1[$i, $i+1] ]); } 
		for (my $i=0; $i<@srt_lis2; $i+=2) { push(@tmp_2, [ @srt_lis2[$i, $i+1] ]); }
		@srt_lis1 = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @tmp_1; 
		@srt_lis2 = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @tmp_2; 
	} elsif ($parm{'sort'} eq '0') {
		my (@tmp_1, @tmp_2); 
		for (my $i=0; $i<@srt_lis1; $i+=2) { push(@tmp_1, [ @srt_lis1[$i, $i+1] ]); } 
		for (my $i=0; $i<@srt_lis2; $i+=2) { push(@tmp_2, [ @srt_lis2[$i, $i+1] ]); }
		@srt_lis1 = @tmp_1; 
		@srt_lis2 = @tmp_2; 
	} else {
		&stopErr("[Err] Unknown 'sort' value [$parm{'sort'}]\n"); 
	}
	
	if ( $parm{'compare'} eq 'same' ) {
		$#srt_lis1 == $#srt_lis2 or return 0; 
		my $str_lis1 = join("\t", map { @$_ } @srt_lis1); 
		my $str_lis2 = join("\t", map { @$_ } @srt_lis2); 
		return ( ($str_lis1 eq $str_lis2) ? 1 : 0 ); 
	} elsif ( $parm{'compare'} eq 'ovl' or $parm{'compare'} eq 'nonovl' ) { 
		my ( $ovlLen, $ovlCnt ) = 0; 
		my @ovlLoc; 
		@srt_lis1 = map { [ sort { $a<=>$b } @$_ ] } @srt_lis1; 
		@srt_lis2 = map { [ sort { $a<=>$b } @$_ ] } @srt_lis2; 
		for my $a1 (@srt_lis1) {
			my ($s1, $e1) = @$a1; 
			for my $a2 (@srt_lis2) {
				my ($s2, $e2) = @$a2; 
				my ($ovl_len, $ovl_loc) = &mathSunhh::ovl_region($s1, $e1, $s2, $e2); 
				$ovl_len > 0 and do { $ovlCnt ++; push( @ovlLoc, [@$ovl_loc] );  }; 
			}
		}
		# Merge @ovlLoc blocks
		@ovlLoc = @{ &mathSunhh::mergeLocBlk(\@ovlLoc) }; 
		
		if ( $parm{'compare'} eq 'ovl' ) {
			for my $a1 (@ovlLoc) {
				$ovlLen += ( $a1->[1] - $a1->[0] + 1 ); 
			}
			return( $ovlLen, $ovlCnt, \@ovlLoc ); 
		} elsif ( $parm{'compare'} eq 'nonovl' ) {
			# Search for non-overlap regions. 
			my @spec1; 
			for my $a1 (@srt_lis1) {
				my ($s1, $e1) = @$a1; 
				for my $a3 (@ovlLoc) {
					my ($s3, $e3) = @$a3; 
					my ( $ovl_o1, $ovl_o2 ) = &mathSunhh::ovl_region($s1, $e1, $s3, $e3); 
					if ( $ovl_o1 > 0 ) {
						$s1 <= $ovl_o2->[0]-1 and push(@spec1, [$s1, $ovl_o2->[0]-1]); 
						$e1 >= $ovl_o2->[1]+1 and push(@spec1, [$ovl_o2->[1]+1, $e1]); 
					} else {
						push(@spec1, [ $s1, $e1 ]); 
					}
					$ovl_o1 > 0 and next; 
				}
			}#End for my $a1 : @spec1
			my @spec2; 
			for my $a1 (@srt_lis2) {
				my ($s1, $e1) = @$a1; 
				for my $a3 (@ovlLoc) {
					my ($s3, $e3) = @$a3; 
					my ( $ovl_o1, $ovl_o2 ) = &mathSunhh::ovl_region($s1, $e1, $s3, $e3); 
					if ( $ovl_o1 > 0 ) {
						$s1 <= $ovl_o2->[0]-1 and push(@spec2, [$s1, $ovl_o2->[0]-1]); 
						$e1 >= $ovl_o2->[1]+1 and push(@spec2, [$ovl_o2->[1]+1, $e1]); 
					} else {
						push(@spec2, [ $s1, $e1 ]); 
					}
					$ovl_o1 > 0 and next; 
				}
			}#End for my $a1 : @spec2
			return ( [@spec1], [@spec2] ); 
		} else {
			&stopErr("[Err] Why here! 'compare'=$parm{'compare'}\n"); 
		}
		return ($ovlLen, $ovlCnt, \@ovlLoc); 
	} else {
		&stopErr("[Err] Unknown 'compare'=[$parm{'compare'}]\n"); 
	}
}# compare_number_list()

=head2 mergeLocBlk( [ [s1,e1], [s2,e2], ... ], [[s21,e21], [s22, e22], ... ] )

Required   : At least one array reference like [ [$s1,$e1] ]; 

Function   : Merge [s1,e1],[s2,e2],[s21,e21],[e22,e22] ... together into non-overlapping blocks. 

 Can be used as &mathSunhh::mergeLocBlk( [ [$s11,$e11], [$s12,$e12] ], [[$s21, $e21]], 'dist2join'=>1  )
 
 Here 'dist2join' is the number distance between two blocks that can be joined. Default is 0, means [4,6] and [7,10] will be two blocks. 
 If set 'dist2join'=>2, [4,6] and [8,10] will be joined into [4,10] ; 

Return     : [[ss1, ee1], [ss2, ee2], ...]

=cut 
sub mergeLocBlk {
	my @back_blk; 
	my @srt_blk; 
	my @para_array; 
	for my $blk_arr (@_) {
		ref($blk_arr) eq '' and do { push(@para_array, $blk_arr); next; }; 
		push(@srt_blk, @$blk_arr); 
	}
	my %parm = &mathSunhh::_setHashFromArr(@para_array); 
	$parm{'dist2join'} //= 0; 
	@srt_blk = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } map { [ sort { $a <=> $b } @$_ ] } @srt_blk; 
	for my $a1 (@srt_blk) {
		my ($s, $e) = @$a1; 
		( defined $s and defined $e ) or next; 
		if ( scalar(@back_blk) > 0 ) {
			if ( $back_blk[-1][1] >= $s - $parm{'dist2join'} ) {
				$e > $back_blk[-1][1] and $back_blk[-1][1] = $e; 
			} else {
				push(@back_blk, [$s, $e]); 
			}
		} else {
			@back_blk = ([$s, $e]); 
		}
	}
	return \@back_blk; 
}# mergeLocBlk() 

=head2 switch_position('qry2ref'=>\%refPos, 'qryID'=>$queryID, 'qryPos'=>$queryPos, 'qryStr'=>'+/-')

Input       : 

 %refPos   : {qryID} => sortedByQrySE [ [qryS_1, qryE_1, refID_1, refS_1, refE_1, refStr_1(+/-)], [qryS_2, qryE_2, refID_2, refS_2, refE_2, refStr_2], ... ]
 
             Requires : qryS <= qryE && refS <= refE 
                        qryS_1 <= qryS_2 ... 

Function    : 

 Get position on ref for qryID_qryPos

Return      : 

 ([refID_1, refPos_1, refStr_1], [refID_2, refPos_2, refStr_2])

=cut
sub switch_position {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	for (qw/qry2ref qryID qryPos/) {
		defined $parm{$_} or &stopErr("[Err] '$_' not defined.\n"); 
	}
	$parm{'qryStr'} //= '+'; 
	my $qry_aref = $parm{'qry2ref'}{ $parm{'qryID'} } // []; 
	my @back; 
	for my $ta (@$qry_aref) {
		$ta->[1] < $parm{'qryPos'} and next; 
		$ta->[0] > $parm{'qryPos'} and last; 
		my $refID = $ta->[2]; 
		my ($refPos, $refStr); 
		$refStr = $parm{'qryStr'}; 
		if ( $ta->[5] eq '+' ) {
			$refPos = $ta->[3] + ($parm{'qryPos'}-$ta->[0]); 
		} elsif ( $ta->[5] eq '-' ) {
			$refPos = $ta->[4] - ($parm{'qryPos'}-$ta->[0]); 
			$refStr =~ tr!+-!-+!; 
		} else {
			&stopErr("[Err] Unknown refStr($ta->[5]) in ID=$parm{'qryID'} [@$ta]\n"); 
		}
		push(@back, [$refID, $refPos, $refStr]); 
	}
	@back == 0 and @back = ([]); 
	return (@back); 
}# sub switch_position () 

=head2 transfer_position( 'from_ref2qry'=>\%agp_scf2ctg, 'to_qry2ref'=>\%agp_ctg2scf, 'fromLoc'=>[$from_scfID, $from_scfPos, $from_scfStrand], 'skipError'=>0 )

Input       : 

  %agp_scf2ctg : Normally it is reversed output of &fileSunhh::load_agpFile(); 
  %agp_ctg2scf : It is directly output of &fileSunhh::load_agpFile();
  Only one-to-one relationship is allowed. 

Return      : return( $to_scfID, $to_scfPos, $to_scfStrand ); 

=cut
sub transfer_position {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	for (qw/from_ref2qry to_qry2ref fromLoc/) {
		defined $parm{$_} or &stopErr("[Err] '$_' not defined in transfer_position().\n"); 
	}
	defined $parm{'fromLoc'}[0] or &stopErr("[Err] No from_scfID defined.\n"); 
	$parm{'fromLoc'}[1] //= 1; 
	$parm{'fromLoc'}[2] //= '+'; 
        $parm{'skipError'}  //= 0;

	my @old_ctgInf = &mathSunhh::switch_position( 'qry2ref' => $parm{'from_ref2qry'} , 'qryID' => $parm{'fromLoc'}[0] , 'qryPos' => $parm{'fromLoc'}[1] , 'qryStr' => $parm{'fromLoc'}[2] ); 
	if (@old_ctgInf == 1 and @{$old_ctgInf[0]} == 0) {
		# &tsmsg("[Wrn] No defined old location for [@{$parm{'fromLoc'}}] [@{$old_ctgInf[0]}]\n"); 
		$old_ctgInf[0] = $parm{'fromLoc'}; 
	} elsif ( @old_ctgInf == 1 and defined $old_ctgInf[0][0] ) {
		; 
        } elsif ( $parm{'skipError'} == 1 ) {
                return('', '', '');
	} else { 
		&stopErr("[Err] Bad result for [@{$parm{'fromLoc'}}] [@old_ctgInf] [@{$old_ctgInf[0]}]\n"); 
	}
	my @new_scfInf = &mathSunhh::switch_position( 'qry2ref' => $parm{'to_qry2ref'} , 'qryID' => $old_ctgInf[0][0] , 'qryPos' => $old_ctgInf[0][1] , 'qryStr' => $old_ctgInf[0][2] );
	if (@new_scfInf == 1 and @{$new_scfInf[0]} == 0) {
		# &tsmsg("[Wrn] No defined new location for [@{$old_ctgInf[0]}]\n"); 
		$new_scfInf[0] = $old_ctgInf[0]; 
	} elsif ( @new_scfInf == 1 ) {
		; 
        } elsif ( $parm{'skipError'} == 1 ) {
                return('', '', '');
	} else { 
		&stopErr("[Err] Bad result for ctg [@{$old_ctgInf[0]}] [@new_scfInf]\n"); 
	}
	return( @{$new_scfInf[0]} ); 
}# transfer_position () 

=head1 sep_loc2_by_loc1_multiLoc2( \@loci_1, \@loci_2, [@colN_from_loci_1] )

Required   : 

 \@loci_1 should be sorted by each start position! 
 \@loci_1 and \@loci_2 in the same format: [ [sec_1_S, sec_1_E], [sec_2_S, sec_2_E], ...  ]

Function   : 

 Separate @loci_2 regions by overlapping @loci_1 or not. 
 And return regions contained in @loci_2; 
 \@loci_2 will be merged into non-overlapping segments. 

Return     : (\@not_overlappedSE, \@overlappedSE)

 Format of @not_overlappedSE

=cut
sub sep_loc2_by_loc1_multiLoc2 {
	my ( $aref_1, $aref_2, $aref_colN ) = @_; 
	my ( @specSE, @ovlpSE, $idx_start ); 
	my $loci_2 = &mathSunhh::mergeLocBlk('', $aref_2, 'dist2join' => 1); 
	for my $ar2 (@$loci_2) {
		my ( $t_specSE, $t_ovlpSE, $t_idx_start ) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( $aref_1, $ar2, $aref_colN ); 
		$idx_start //= $t_idx_start; 
		$idx_start > $t_idx_start and $idx_start = $t_idx_start; 
		push(@specSE, @$t_specSE); 
		push(@ovlpSE, @$t_ovlpSE); 
	}
	
	return( \@specSE, \@ovlpSE, $idx_start ); 
}# sep_loc2_by_loc1_multiLoc2() 

=head1 sep_loc2_by_loc1_singleLoc2( \@loci_1, [$s2, $e2], [@colN_from_loci_1] ) 

Input      : 

 \@loci_1 should be sorted by each start position! 
 \@loci_1 in the format: [ [sec_1_S, sec_1_E], [sec_2_S, sec_2_E], ...  ]
 [@colN_from_loci_1] : Default [], means no columns to add. 
   If given, like [0,3], the 0-th and 3-rd columns from each $loci_1[$x] element will be added to \@overlappedSE; 

Required   : \@loci_1 , [$s2, $e2]

Function   : 

 Separate [$s2, $e2] region by overlapping @loci_1 or not. 
 And return regions contained in [$s2, $e2]

Return     : (\@not_overlappedSE, \@overlappedSE, $idx_start)

 Format of @not_overlappedSE : ( [ $s_1, $e_1 ], [ $s_2, $e_2 ], ... )
 Format of @overlappedSE : ( [ $s_1, $e_1, ... added ], [ $s_2, $e_2, ... added ], ... )
 $idx_start : The elements before this index_element is not larger than any of \@not_overlappedSE|\@overlappedSE . 

=cut
sub sep_loc2_by_loc1_singleLoc2 {
	my ($aref_1, $aref_2, $aref_colN) = @_; 
	my (@specSE, @ovlpSE); 
	$aref_colN //= []; 
	
	my $cur_p = $aref_2->[0]; 
	my $idx_cnt = -1; 
	my $idx_start; #  
	for my $ar1 (@$aref_1) {
		# ar1 -> ( loci_1_start, loci_1_end )
		$idx_cnt ++; 
		$ar1->[1] < $cur_p and next;       # The end of loci_1 (end_1)is smaller than current positoin; Go to next. 
		$ar1->[0] > $aref_2->[1] and last; # The beginning of loci_1 (start_1) is larger than loci_2, so we don't need to compare the followings. 
		if ( $ar1->[0] <= $cur_p) {
			# The start_1 is smaller than current position, but end_1 is no less than current position; 
			# The loci_1 is spanning current position. 
			$idx_start //= $idx_cnt; 
			if ( $ar1->[1] >= $cur_p ) {
				if ( $ar1->[1] <= $aref_2->[1] ) {
					push(@ovlpSE, [ $cur_p, $ar1->[1], @{$ar1}[@$aref_colN] ]); 
				} else {
					push(@ovlpSE, [ $cur_p, $aref_2->[1], @{$ar1}[@$aref_colN] ]); 
				}
			}
			$cur_p = $ar1->[1] + 1; 
		} elsif ( $ar1->[0] > $cur_p ) {
			# The loci_1 is overlapping loci_2, but larger than current position
			$idx_start //= $idx_cnt; 
			push( @specSE, [$cur_p, $ar1->[0]-1 ] ); 
			$cur_p = $ar1->[1] + 1; 
			if ( $ar1->[1] <= $aref_2->[1] ) {
				push( @ovlpSE, [ $ar1->[0], $ar1->[1], @{$ar1}[@$aref_colN] ] ); 
			} else {
				push( @ovlpSE, [ $ar1->[0], $aref_2->[1], @{$ar1}[@$aref_colN] ] ); 
			}
		} else {
			# This is impossible! 
			&stopErr("[Err] Wrongly reach this line in sep_loc2_by_loc1_singleLoc2()\n"); 
		}
	}# for my $ar1 (@$aref_1)
	
	if ( $cur_p <= $aref_2->[1] ) {
		# After checking loci_1 regions, there is also loci_2 region left. 
		push( @specSE, [ $cur_p, $aref_2->[1] ] ); 
		$idx_cnt < 0 and $idx_cnt = 0; # In this case, the loci_1 is empty. 
		$idx_start //= $idx_cnt; 
	}
	
	return (\@specSE, \@ovlpSE, $idx_start); 
}# sep_loc2_by_loc1_singleLoc2() 

=head1 map_loc_byReference( \@reference_1to2, \@LocSE_1, '+/-' )

Input      : 

 \@reference_1to2 : Format as [ $start1, $end1, $start2, $end2 ], in which loc_1 [$start1,$end1] is mapped to loc_2 [$start2, $end2]; 
 \@LocSE_1        : Format as [ $targetS_1, $targetE_1 ] 
 '+/-'            : The strand (direction) how loc_1 is aligned to loc_2. If '-', $start1 is mapped to $end2, and $end1 is mapped to $start2; 

Function   : 

 According to \@reference_1to2, map loc_1 positions [ $targetS_1, $targetE_1 ] to loc_2 positions, and return them. 

Return     : ( $targetS_2, $targetE_2 )

=cut
sub map_loc_byReference {
	my ( $aref_1, $aref_2, $str ) = @_; 
	$str //= '+'; 
	my ( $targetS_2, $targetE_2 ); 
	
	# $aref_1 is in format [ $start1, $end1, $start2, $end2 ]
	# $aref_2 is in format [ $targetS_1, $targetE_1 ]
	my $ratio = ( $aref_1->[1]-$aref_1->[0] == 0 ) ? 0 : ($aref_1->[3]-$aref_1->[2])/($aref_1->[1]-$aref_1->[0]) ; 
	## Here $targetS_2 should be always no larger than $targetE_2, in case $aref_2->[1] >= $aref_2->[0]; 
	$targetS_2 = ( $str eq '-' ) 
	  ? $aref_1->[3] - int( ($aref_2->[1]-$aref_1->[0]) * $ratio)
	  : $aref_1->[2] + int( ($aref_2->[0]-$aref_1->[0]) * $ratio )
	; 
	$targetE_2 = ( $str eq '-' )
	  ? $aref_1->[3] - int( ($aref_2->[0]-$aref_1->[0]) * $ratio ) 
	  : $aref_1->[2] + int( ($aref_2->[1]-$aref_1->[0]) * $ratio )
	; 
	
	return ($targetS_2, $targetE_2); 
}# map_loc_byReference() 

=head1 insert_SEloc_toSEArray( \@toSEArray, \@SEloc, $idx_start )

Input      : 

 \@toSEArray   : Sorted. [ [$s1,$e1], [$s2, $e2], ... ]
 \@SEloc       : [ $s, $e ]
 $idx_start    : The elements in @toSEArray before index($idx_start) is not smaller than $s. Default 0. 

Required   : \@toSEArray, \@SEloc

Function   : Insert \@SEloc into \@toSEArray by low-to-high order. 

Return     : No return. The array \@toSEArray will be changed. 

=cut
sub insert_SEloc_toSEArray {
	my ($aref_1, $aref_2, $idx_curr) = @_; 
	$idx_curr //= 0; 
	my $is_found = 0; 
	for (; $idx_curr < @$aref_1; $idx_curr++) {
		my $ar1 = $aref_1->[$idx_curr]; 
		# $ar1 -> [ start, end ]
		$ar1->[0] < $aref_2->[0] and next; 
		$ar1->[0] == $aref_2->[0] and $ar1->[1] <= $aref_2->[1] and next; 
		$is_found = 1; 
		last; 
	}
	if ( $is_found == 0 or $idx_curr > $#$aref_1 ) {
		push(@$aref_1, $aref_2); 
	} else {
		splice(@$aref_1, $idx_curr, 0, $aref_2); 
	}
	return; 
}# insert_SEloc_toSEArray() 

1; 
