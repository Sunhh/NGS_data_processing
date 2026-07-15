package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 


=head1 index_SEloc( [ [$loc_1_S, $loc_1_E], [$loc_2_S, $loc_2_E], ... ], {'binSize' => 1000 , 'maxP' => 11 } )

Return        : (\%db_loc)

Usage case 1  : 

    my @loc = ( [1,100], [2 , 384745], [10,21], [50,60], [384745, 984672], [20,55], [1, 9999], [4,11], [10,20], [10,21], [10,19]); 
    my %locDb_h = %{ &mathSunhh::index_SEloc( \@loc, { 'binSize' => 100 , 'maxP' => 11 } ) }; 
    print "intact raw data Start:\n"; 
    for (my $i=0; $i<@loc; $i++) {
      print join("\t", "$loc[$i][0]-$loc[$i][1]", "$locDb_h{'realLoc'}[$i][0]-$locDb_h{'realLoc'}[$i][1]")."\n"; 
    }
    print "intact raw data End:\n"; 


    ## Check if it is well sorted : 
    print "newly sorted Start :\n"; 
    my @sorted_realLoc = @{ &mathSunhh::ret_sortedRealLoc( \%locDb_h ) }; 
    for my $ar1 (@sorted_realLoc) {
      print join("\t", @$ar1)."\n"; 
    }
    print "newly sorted End :\n"; 

    ## Find overlapping realLocs with given SEloc
    my @seLoc = (21, 30); 
    @ARGV >= 2 and @seLoc = @ARGV[0,1]; 
    print "\nRegion to overlap to: $seLoc[0] - $seLoc[1]\n\n"; 
    my @loc_realIdx = @{ &mathSunhh::_map_loc_to_realIdx( \%locDb_h, \@seLoc ) }; 
    my @loc_realLoc = @{$locDb_h{'realLoc'}}[@loc_realIdx]; 
    for (my $i=0; $i<@loc_realIdx; $i++) {
      print join("\t", "realIdx:$loc_realIdx[$i]", "realLoc:_$loc_realLoc[$i][0]_-_$loc_realLoc[$i][1]")."\n";
    }



Description   : 

    This is the main function I want to use in RNAseq counting. 
    $loc_N_S should be no larger than (<=) $loc_N_E ; 

 By self()
    Add ->{'binSize'}             => 1000 ; 
    Add ->{'maxP'}                => 11 ; 
    Add ->{'realLoc'}             => [ [$rawLoc_1_S, $rawLoc_1_E], [$rawLoc_2_S, $rawLoc_2_E], ... ]
    Add ->{'minV'}                => $min_of_locS_locE_value; 
    Add ->{'maxV'}                => $max_of_locS_locE_value; 
 By _binRealLoc2BinLoc() 
    Add ->{'rawIdx2realIdx'}{$rawIdx}    -> [$realIdx_1, $realIdx_2, ...]
          # This relationship is the most important to transform 'rawLoc' to 'realLoc' information. 
    Add ->{'rawLoc'}              => [ [$rawLoc_1_S, $rawLoc_1_E], [$rawLoc_2_S, $rawLoc_2_E], ... ]
          # Here 'rawLoc' stands for 'binLoc'; 
    Add ->{'rawNum'}              => [ $rawLoc_1_S, $rawLoc_2_S, $rawLoc_3_S, ... ]
          # Here 'rawNum' is related to 'rawLoc' instead of 'realLoc'
 By _fill_rawN2rawIdx()
    Add ->{'rawN2rawIdx'}{$rawLoc_S} => [ $rawIdx1, $rawIdx2, ... ]
                         # Here $rawLoc_S === $rawNum
                         # This value will be replaced by the following _fill_struct_SEloc() function. 
 By _fill_struct_SEloc()
    Add ->{'struct'}           => [ structure_database ]
                                  # final value of 'struct' is ->{'rawN2rawIdx'}{$rawN} in a corrected order. 
                                  # and here $rawN is the start of each location. 
                                  # Formatted as [ $rawIdx1, $rawIdx2, ... ]
                                  # Note that the 'raw' here related to 'rawLoc' instead of 'realLoc'. 
    Replace ->{'rawN2rawIdx'}  => $final_value_of_"$ah->{'struct'}", which is now in a corrected order. 
 By _fill_sorted_rawIdx()
    Add ->{'sorted_rawIdx'} => [ [$rawIdx3, $rawIdx8, ...], [$rawIdx2], ... ];
          # Here $rawIdx3_rawNumber is equal to $rawIdx8_rawNumber, and smaller than $rawIdx2_rawNumber;
          # The array_reference in values is linked to ->{'rawN2rawIdx'}{$rawN} for easy changes;
    Add ->{'rawN2srtIdx'}{$rawN} => $sorted_idx_in_'sorted_rawIdx_array';

=cut
sub index_SEloc {
	my ($arLoc, $ph) = @_; 
	my %back; 
	$ph->{'binSize'} //= 1000; $back{'binSize'} = $ph->{'binSize'}; 
	$ph->{'maxP'}    //= 11;   $back{'maxP'}    = $ph->{'maxP'}; 
	@{$back{'realLoc'}} = map { ($_->[0] > $_->[1]) ? [ @{$_}[1,0] ] : [ @{$_}[0,1] ] ; } @$arLoc; 
	($back{'minV'}, $back{'maxV'}) = &mathSunhh::minmax( [ map { ($_->[0], $_->[1]) } @{$back{'realLoc'}} ] ); 
	&mathSunhh::_binRealLoc2BinLoc( \%back ); 
	&mathSunhh::_fill_rawN2rawIdx( \%back ); 
	&mathSunhh::_fill_struct_SEloc( \%back ); 
	&mathSunhh::_fill_sorted_rawIdx( \%back ); 
	
	return(\%back); 
}# index_SEloc ()

=head1 _map_loc_to_realIdx ( \%db_loc, [$chkLoc_S, $chkLoc_E] ) 

Return       : ( [ $overlap_realLocIdx_1, $overlap_realLocIdx_2, ... ] )

Description  : $overlap_realLocIdx can be used in $db_loc{'realLoc'}[$overlap_realLocIdx] (Value=[$realS, $realE]); 

=cut
sub _map_loc_to_realIdx {
	my ($locDB_hr, $locSE_ar) = @_; 
	my @back; 
	my @loc_rawIdx = @{ &mathSunhh::_map_loc_to_rawIdx( $locDB_hr , $locSE_ar ) }; 
	my ($tS,$tE) = @$locSE_ar; 
	my @realIdx_arr = @{ &mathSunhh::_locDb_rawIdx2realIdx( $locDB_hr, \@loc_rawIdx, 1 ) }; 
#	# I don't want to make it too too too complex, so here I just compare loci from the beginnning. Or esle I should try to use _closest_number(); 
#	# This loci_array should be from small to big. (increasing)
	for my $realIdx ( @realIdx_arr ) {
		my ($curS, $curE) = @{ $locDB_hr->{'realLoc'}[$realIdx] }; 
		$curS > $tE and last; 
		$curE < $tS and next; 
		push(@back, $realIdx); 
	}

	return(\@back); 
} # _map_loc_to_realIdx () 

=head1 ret_sortedRealLoc (\%db_loc) 

Return     : ( [ [loc1_S, loc1_E], [loc2_S, loc2_E], ... ] )

=cut
sub ret_sortedRealLoc {
	my ($ah) = @_; 
	my @back; 
	my @sorted_rawIdx = map { @$_ } @{$ah->{'sorted_rawIdx'}}; 
	my @sorted_realIdx = @{ &mathSunhh::_locDb_rawIdx2realIdx( $ah, \@sorted_rawIdx, 1 ) }; 
	@back = map { [@$_] } @{$ah->{'realLoc'}}[ @sorted_realIdx ]; 

	return(\@back); 
}# ret_sortedRealLoc () 

=head1 index_numbers( \@raw_numbers, { 'maxP' => 11 } )

Return       : (\%db_num)

Function     : 

 Index the numbers, producing keys in \%db_num, including : {'maxP|rawNum|rawN2rawIdx|struct|sorted_rawIdx|rawN2srtIdx'}
    ->{'maxP'}          =>  11 , or someother number (>0) showing the depth of structure and the depth (length)  of encoded number path (for &mathSunhh::_Encode_deciN())
    ->{'rawNum'}        =>  [@raw_numbers]
    ->{'rawN2rawIdx'}   =>  { $rawNum => [ $rawIdx1, $rawIdx2, ... ] } according to ->{'rawNum'}
    ->{'struct'}        =>  An array of arrays with depth (rows) of ->{'mapP'} and width (cols) of [0-9] or the maximum of the first(0) row. 
    ->{'sorted_rawIdx'} =>  [ ->{'rawN2rawIdx'}{$rawNum_1}, {'rawN2rawIdx'}{$rawNum_2}, ... ] 
                            === [ [$rawIdx3, $rawIdx8, ...], [$rawIdx2], ... ], 
                               # Here, $rawIdx3_rawNumber ( === ->{'rawNum'}[$rawIdx3] ) is equal to $rawIdx8_rawNumber, and smaller than $rawIdx2_rawNumber;
    ->{'rawN2srtIdx'}   =>  { $rawNum => $index_of_'sorted_rawIdx'_array }

 Usage case 1 : 
    my @numbers = ( 1.. 100, 40 .. 50 ,105, 103, 105, 150, 1, 999, 3000, 888, 444, 999, 3000, 999, 99999999999999, 111, 99999999999399 ); 
    my %numDb_h = %{ &mathSunhh::index_numbers( \@numbers, { 'maxP' => 11 } ) }; 

    print "input\@=@numbers\n"; 
    print "stored=@{$numDb_h{'rawNum'}}\n"; 
    ## Get sorted numbers : 
    print "sorted=" . join(" ", @{ &mathSunhh::ret_sortedNum( \%numDb_h )}) . "\n"; 
    ## Get sorted raw indices : 
    my @sorted_rawIdx = map { @$_ } @{$numDb_h{'sorted_rawIdx'}}; 
    for my $rawIdx ( @sorted_rawIdx ) {
      print join("\t", "rawIdx:$rawIdx", "rawNum:$numDb_h{'rawNum'}[$rawIdx]")."\n"; 
    }

    ## Locate a number in database : 
    my $number_to_check = 99999999999300; 
    print "searching=$number_to_check\n"; 
    my $srtIdx = &mathSunhh::_map_number_to_srtIdx( \%numDb_h, $number_to_check ); 
    print join("\t", "srchN:$number_to_check", "srtIdx:$srtIdx", join(":", "rawNumbers", @{$numDb_h{'rawNum'}}[ @{$numDb_h{'sorted_rawIdx'}[$srtIdx]} ]))."\n";
    # Printing : srchN:10000     srtIdx:107      rawNumbers:3000:3000

=cut
sub index_numbers {
	my ($arN, $pH) = @_; 
	$pH->{'maxP'} //= 11; 
	my %back; 
	$back{'maxP'} = $pH->{'maxP'}; 
	### Add ->{'rawNum'}        => [ $rawN_1, $rawN_2, ... ]
	@{$back{'rawNum'}} = @$arN; 
	&mathSunhh::_fill_rawN2rawIdx( \%back ); 
	&mathSunhh::_fill_struct( \%back ); 
	&mathSunhh::_fill_sorted_rawIdx( \%back ); 

	return(\%back); 
}# sub index_numbers() 


=head1 _sideClose_pathInStruct( \@struct, \@path, 'both|left|right' )

Return      : ( \@path_of_closest )

When using 'both' sides, I will return the smaller one if the two sides are equally distant. 

=cut
sub _sideClose_pathInStruct {
	my ( $ar_struct, $ar_path_in, $side ) = @_; 
	$side //= 'both'; 
	$side =~ m/^(left|right|both)$/ or &stopErr("[Err] Bad side [$side] parameter in _sideClose_pathInStruct()\n"); 
	my $path_depth = $#$ar_path_in + 1; 

	my $tmp_ar = $ar_struct; 
	for (my $i=0; $i<$path_depth; $i++) {
		if ( defined $tmp_ar->[ $ar_path_in->[$i] ] ) {
			$tmp_ar = $tmp_ar->[ $ar_path_in->[$i] ]; 
			next; 
		}
		if ( $side eq 'left' ) {
			# S1: Check if the current level ($i) has a idxV < $ar_path_in->[$i]; 
			for (my $j=$ar_path_in->[$i]-1; $j>=0; $j--) {
				defined $tmp_ar->[$j] or next; 
				# If has such an element 
				return( &mathSunhh::_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $i-1], $j ], 'max', $path_depth ) ); 
			}
			# S2: There is no idxV smaller than input in current level ($i), so find the biggest in the most recent uppper level. 
			#     I don't want to combine this and the S1 step, because I want to use &mathSunhh::_defined_pathInStruct() to confirm there is no problem. 
			for (my $lvl=$i-1; $lvl >= 0; $lvl--) {
				my $lvl_struct = &mathSunhh::_defined_pathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl] ] ) or &stopErr("[Err] Weired get here! 'left'\n"); 
				my $lvl_idx    = $ar_path_in->[$lvl]; 
				# Check if this level ($lvl) has an idxV < $ar_path_in->[$lvl]; 
				for ( my $j=$lvl_idx-1; $j>=0; $j-- ) {
					defined $lvl_struct->[$j] or next; 
					return( &mathSunhh::_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl], $j ], 'max', $path_depth ) ); 
				}
			}
			# S3: There is no smaller path found in total dataset, so I should use the smallest in database. 
			return( &mathSunhh::_get_pathInStruct( $ar_struct, 'min', $path_depth ) ); 
		} elsif ( $side eq 'right' ) {
			# S1: Check if the current level ($i) has a idxV > $ar_path_in->[$i]; 
			for (my $j=$ar_path_in->[$i]+1; $j<=@$tmp_ar; $j++) {
				defined $tmp_ar->[$j] or next; 
				# If has such an element
				return( &mathSunhh::_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $i-1], $j ], 'min', $path_depth ) ); 
			}
			# S2: Continue to check from the most recent upper level ($i-1). 
			for (my $lvl=$i-1; $lvl >= 0; $lvl--) {
				my $lvl_struct = &mathSunhh::_defined_pathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl] ] ) or &stopErr("[Err] Wiered get here! 'right'\n"); 
				my $lvl_idx    = $ar_path_in->[$lvl]; 
				# Check if this level contains an idxV > $ar_path_in->[$lvl]; 
				for ( my $j=$lvl_idx+1; $j<=@$lvl_struct; $j++ ) {
					defined $lvl_struct->[$j] or next; 
					return( &mathSunhh::_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl], $j ], 'min', $path_depth ) ); 
				}
			}
			# S3: There is no smaller path found in total dataset, so I should use the biggest in database. 
			return( &mathSunhh::_get_pathInStruct( $ar_struct, 'max', $path_depth ) ); 
		} elsif ( $side eq 'both' ) {
			# S1: Find the left most path 
			my $path_left = &mathSunhh::_sideClose_pathInStruct( $ar_struct, $ar_path_in, 'left' ); 
			my $dist_left = abs( &mathSunhh::_distance_paths( $ar_path_in, $path_left ) ); 
			# S2: Find the right most path 
			my $path_right = &mathSunhh::_sideClose_pathInStruct( $ar_struct, $ar_path_in, 'right' ); 
			my $dist_right = abs( &mathSunhh::_distance_paths( $path_right, $ar_path_in ) ); 
			# S3: Compare distances between chk-left and right-chk, and find the smaller one. 
			if ( $dist_left <= $dist_right ) {
				return( $path_left ); 
			} else {
				return( $path_right ); 
			}
		} else {
			&stopErr("[Err] Bad side information [$side]\n"); 
		}
	}

	# Arrive here only when the $ar_path_in exists. 
	return([@$ar_path_in]); 
}# _sideClose_pathInStruct ()

=head1 _distance_paths( $path_1, $path_2 )

Return        : ( $diff_value_by_decode_num )

 num($path_1) - num($path_2)

=cut
sub _distance_paths{
	my ($p1, $p2) = @_; 
	my $back = &mathSunhh::_Decode_deciN($p1) - &mathSunhh::_Decode_deciN($p2); 
	return($back); 
}# _distance_paths () 

=head1 _defined_pathInStruct( \@struct, \@path ) 

Return        : ('0|$ar_sub_struct') # 0 - undefined. 

=cut
sub _defined_pathInStruct {
	my ($ar_struct, $path) = @_; 
	my $tmp_ar = $ar_struct; 
	for (my $i=0; $i<@$path; $i++) {
		$path->[$i] < 0 and return 0; 
		defined $tmp_ar->[$path->[$i]] or return 0; 
		$tmp_ar = $tmp_ar->[$path->[$i]]; 
	}
	return $tmp_ar; 
}# _defined_pathInStruct () 

=head1 _lsAll_inStruct( \@struct, $maxP )

Return        : ( [ [$struct_path, $struct_final_value], [$struct_path, $struct_final_value], ...  ] )

 From small (first) to big (last) order. 

=cut
sub _lsAll_inStruct {
	my ($ar_struct, $maxP) = @_; 
	my @back; 
	$maxP <= 0 and &stopErr("[Err] maxP must be bigger than 0.\n"); 
	if ($maxP == 1) {
		for (my $i=0; $i<@$ar_struct; $i++) {
			defined $ar_struct->[$i] or next; 
			push(@back, [ [$i], $ar_struct->[$i] ]); 
		}
		return(\@back); 
	} elsif ( $maxP > 1 ) {
		for (my $i=0; $i<@$ar_struct; $i++) {
			defined $ar_struct->[$i] or next; 
			my $sub_back = &mathSunhh::_lsAll_inStruct( $ar_struct->[$i], $maxP-1 ); 
			for my $ar1 ( @$sub_back ) {
				push(@back, [ [ $i, @{$ar1->[0]} ], $ar1->[1] ]); 
			}
		}
		return(\@back); 
	} else {
		&stopErr("[Err] Bad here in _lsAll_inStruct()\n"); 
	}
	return; 
}# _lsAll_inStruct () 

=head1 ret_sortedNum (\%db_number) 

Return       : ( [smallest_number1, smallest_number2, ...] )

=cut
sub ret_sortedNum {
	my ($ah) = @_; 
	return( [ map { @{ $ah->{'rawNum'} }[@$_] } @{$ah->{'sorted_rawIdx'}} ] ); 
}# ret_sortedNum () 

=head1 _map_number_to_rawNum( \%db_number, $number_to_locate, 'both|left|right' )

Return       : ($rawNumber)

Description  : The third parameter will be 'both' in default. 

=cut
sub _map_number_to_rawNum {
	my ($ah, $number, $side) = @_; 
	$side //= 'both'; 
	my $num_path = &mathSunhh::_Encode_deciN($number); 
	my $closest_path = &mathSunhh::_sideClose_pathInStruct( $ah->{'struct'}, $num_path, $side ); 
	my $closest_rawN = $ah->{'rawNum'}[ &mathSunhh::_get_valueFromStruct( $ah->{'struct'}, $closest_path )->[0] ]; 

	return( $closest_rawN ); 
}# _map_number_to_rawNum ()

=head1 _map_number_to_srtIdx( \%db_number, $number_to_locate, 'both|left|right' )

Return       : ($idx_in_'sorted_rawIdx_array')

Description  : The third parameter will be 'both' in default. 

=cut
sub _map_number_to_srtIdx {
	my ($ah, $number, $side) = @_; 
	$side //= 'both'; 
	my $closest_rawN = &mathSunhh::_map_number_to_rawNum( $ah, $number, $side ); 
	# my $closest_rawN = &mathSunhh::_Decode_deciN($closest_path) ; 
	my $closest_srtIdx = $ah->{'rawN2srtIdx'}{$closest_rawN}; 

	return( $closest_srtIdx ); 
}# _map_number_to_srtIdx ()

=head1 _get_valueFromStruct( \@struct, \@encoded_num )

Return        : ($final_value)

Function      : 
                $final_value should be $num_db_hash{'rawN2rawIdx'}{$rawN}; 

=cut
sub _get_valueFromStruct {
	my ( $ar_1, $ar_2 ) = @_; 
	my $back; 
	my $tmp_ar = $ar_1; 
	for (my $i=0; $i<@$ar_2; $i++) {
		$tmp_ar = $tmp_ar->[$ar_2->[$i]]; 
	}
	$back = $tmp_ar; 

	return($back); 
}# _get_valueFromStruct () 

=head1 _get_subPathInStruct( \@struct, \@root_sub_path, 'min|max', $path_depth ) 

Return        : (\@path_of_minORmax)

=cut
sub _get_subPathInStruct {
	my ( $ar_struct, $root_path, $mm, $path_depth ) = @_; 
	my @back; 
	$mm //= 'min'; 
	$mm =~ m/^(min|max)$/ or &stopErr("[Err] bad mm_type [$mm] in _get_subPathInStruct()\n"); 

	my $sub_struct = &mathSunhh::_defined_pathInStruct( $ar_struct, $root_path ) or &stopErr("[Err] No structure exists\n"); 
	my $sub_path   = &mathSunhh::_get_pathInStruct( $sub_struct, $mm, $path_depth-scalar(@$root_path) ); 
	@back = ( @$root_path, @$sub_path ); 
	scalar(@back) == $path_depth or &stopErr("[Err] Unequal depth of resulting path ($#back+1) and input path depth ($path_depth).\n"); 

	return(\@back); 
}# _get_subPathInStruct () 

=head1 _get_pathInStruct( \@struct, 'min|max', $path_depth )

Return        : (\@path_in_struct)

Function      : Get the min/max path in structure. All three parameters needed. 

=cut
sub _get_pathInStruct {
	my ( $ar_struct, $type, $path_depth ) = @_; 
	$path_depth //= 11; 
	my $tmp_ar = $ar_struct; 
	my @back = (0) x $path_depth; 

	if ( $type eq 'min' or $type eq 'smallest' ) {
#		# This method works in case that we know the structure of the final value . 
#		my $find = 1; 
#		while ( defined $tmp_ar and $find ) {
#			$find = 0; 
#			for (my $i=0; $i<@$tmp_ar; $i++) {
#				defined $tmp_ar->[$i] or next; 
#				$tmp_ar = $tmp_ar->[$i]; 
#				push(@back, $i); 
#				$find = 1; 
#				last; 
#			}
#		}
		for (my $i=0; $i<$path_depth; $i++) {
			for (my $j=0; $j<@$tmp_ar; $j++) {
				defined $tmp_ar->[$j] or next; 
				$tmp_ar = $tmp_ar->[$j]; $back[$i] = $j; 
				last; 
			}
		}
	} elsif ( $type eq 'max' or $type eq 'biggest' ) {
		for (my $i=0; $i<$path_depth; $i++) {
			for (my $j=$#$tmp_ar; $j>=0; $j--) {
				defined $tmp_ar->[$j] or next; 
				$tmp_ar = $tmp_ar->[$j]; $back[$i] = $j; 
				last; 
			}
		}
	} else {
		&stopErr("[Err] Bad type [$type] for _get_pathInStruct()\n"); 
	}

	return(\@back); 
}# _get_pathInStruct()

=head1 _put_num2Struct( \@struct, \@encoded_num, $final_value )

Return        : ()

Function      : 

                Change values in @struct, which will be used in @{$num_db_hash{'struct'}}; 
                @encoded_num comes from &mathSunhh::_Encode_deciN() function, and it is treated as a full path in @struct. 
                $final_value will be used for the terminal element if that element has not been defined. 
                And $final_value should be generated as [ @{$num_db_hash{'rawN2rawIdx'}{$rawN}} ]. 
                Here I directly use $final_value in order to make linkage between this element and $num_db_hash{'rawN2rawIdx'}{$rawN}; 

=cut
sub _put_num2Struct {
	my ( $ar_1, $ar_2, $ar_3 ) = @_; 
	# @$ar_1 is database structure. 
	# @$ar_2 is the path from number (encoded by _Encode_deciN() )
	# @$ar_3 is the array for value set. 
	my $tmp_ar = $ar_1; 
	for (my $i=0; $i<@$ar_2; $i++) {
		$tmp_ar->[$ar_2->[$i]] //= ( ( $i == $#$ar_2 ) ? $ar_3 : [] ); 
		$tmp_ar = $tmp_ar->[$ar_2->[$i]]; 
	}

	return; 
}# _put_num2Struct () 

=head1 _locDb_rawIdx2realIdx( \%db_loc, \@rawIdx_ar, [1|0"if removing duplicated realIdx"] )

Return        : ( [realLocIdx_1, realLocIdx_2, ...] )

Description   : 

 Change 'rawIdx' in database to 'realIdx' numbers; 
 If the third parameter is true (default) (not zero or null), the duplicated realIdx will be removed. 

=cut
sub _locDb_rawIdx2realIdx {
	my ($ah, $rawIdx_ar, $rmDup) = @_; 
	$rmDup //= 1; 
	my @back; 
	my %used_realIdx; 
	for my $realIdx ( map { @{ $ah->{'rawIdx2realIdx'}{$_} } } @$rawIdx_ar ) {
		defined $used_realIdx{$realIdx} and next; 
		$rmDup and $used_realIdx{$realIdx} //= 1; 
		push(@back, $realIdx); 
	}
	return(\@back); 
}# _locDb_rawIdx2realIdx() 







############################################################
#  Inner sub-routines. 
############################################################

=head1 _fill_rawN2rawIdx( \%num_db_hash ) 

Return        : () 

Function      : 

    Add ->{'rawN2rawIdx'}{$rawNum} => [ $rawIdx1, $rawIdx2, ... ] according to ->{'rawNum'}

=cut
sub _fill_rawN2rawIdx {
	my ($ah) = @_; 
	for (my $i=0; $i<@{ $ah->{'rawNum'} }; $i++) {
		push(@{$ah->{'rawN2rawIdx'}{ $ah->{'rawNum'}[$i] }}, $i); 
	}
	return; 
}# _fill_rawN2rawIdx ()

=head1 _fill_struct( \%db_num ) 

Function    : 

 Add ->{'struct'}     => [ structure_database ]
 final value of 'struct' is ->{'rawN2rawIdx'}{$rawN} 

=cut
sub _fill_struct {
	my ( $ah ) = @_; 
	$ah->{'struct'} = []; 
	for my $rawN (keys %{$ah->{'rawN2rawIdx'}}) {
		&mathSunhh::_put_num2Struct( $ah->{'struct'}, &mathSunhh::_Encode_deciN($rawN, $ah->{'maxP'}), $ah->{'rawN2rawIdx'}{$rawN} ); 
	}
	return; 
}# _fill_struct () 

=head1 _fill_struct_SEloc( \%db_loc , { 'maxP' => 11 } )

Return       : ()

Function     : 

    The second added 'maxP' only works when $db_loc{'maxP'} is undefined. 
    Add ->{'struct'}           => [ structure_database ]
                                  # final value of 'struct' is ->{'rawN2rawIdx'}{$rawN} in a corrected order. 
                                  # and here $rawN is the start of each location. 
    Replace ->{'rawN2rawIdx'}  => $final_value_of_"$ah->{'struct'}", which is now in a corrected order. 

=cut
sub _fill_struct_SEloc {
	my ($ah, $pH) = @_; 
	$pH->{'maxP'} //= 11; 
	$ah->{'maxP'} //= $pH->{'maxP'}; 

	# Firstly sorted by locE : 
	my @eLoc = map { $_->[1] } @{$ah->{'rawLoc'}}; 
	my %eLoc_db = %{ &mathSunhh::index_numbers(\@eLoc, { 'maxP' => $ah->{'maxP'} }) }; 
	my %eLoc_info; 
	for my $eLoc_rawIdx_ar ( @{$eLoc_db{'sorted_rawIdx'}} ) {
		for my $eLoc_rawIdx_v (@$eLoc_rawIdx_ar) {
			# Because @eLoc directly comes from @{$ah->{'rawLoc'}}, so $eLoc_rawIdx_v === ah:rawLocIdx; 
			push(@{$eLoc_info{'srt_locSE'}}, $ah->{'rawLoc'}[$eLoc_rawIdx_v]); 
			$eLoc_info{'srtIdx_to_rawIdx'}{ $#{$eLoc_info{'srt_locSE'}} } = $eLoc_rawIdx_v; # This srtIdx is not index of 'sorted_rawIdx'!!! 
		}
	}

	# Secondly sorted by locS : 
	my @sLoc = map { $_->[0] } @{$eLoc_info{'srt_locSE'}}; 
	my %sLoc_db = %{ &mathSunhh::index_numbers(\@sLoc, { 'maxP' => $ah->{'maxP'} } ) }; 
	my %sLoc_info; 
	for my $sLoc_rawIdx_ar ( @{$sLoc_db{'sorted_rawIdx'}} ) {
		for my $sLoc_rawIdx_v (@$sLoc_rawIdx_ar) {
			# Here $sLoc_rawIdx_v === $eLoc_srtIdx ; 
			my $eLoc_srtIdx = $sLoc_rawIdx_v; 
			push(@{$sLoc_info{'srt_locSE'}}, [ @{$eLoc_info{'srt_locSE'}[$eLoc_srtIdx]} ]); 
			my $sLoc_srtIdx = $#{$sLoc_info{'srt_locSE'}}; 
			# my $rawIdx = $eLoc_info{'srtIdx_to_rawIdx'}{$eLoc_srtIdx}; 
			# $sLoc_info{'srtIdx_to_rawIdx'}{$sLoc_srtIdx} = $rawIdx; 
		}
	}
	
	# Now we can fill 'struct' related values : 
	$ah->{'struct'} = [ @{ $sLoc_db{'struct'} } ]; 
	for my $ar_rawIdx_INsLocDb ( map { $_->[1] } @{ &mathSunhh::_lsAll_inStruct( $ah->{'struct'} , $ah->{'maxP'} ) } ) {
		my $sLoc_rawLocS = $sLoc_db{'rawNum'}[ $ar_rawIdx_INsLocDb->[0] ]; 
		@$ar_rawIdx_INsLocDb = @{$eLoc_info{'srtIdx_to_rawIdx'}}{@$ar_rawIdx_INsLocDb}; # Remember $sLoc_rawIdx_v === $eLoc_srtIdx ;
		# We want to renew 'rawN2rawIdx' information, because we want to sort loci also with Eloc. 
		# We have to do this because the current ->{'struct'} is constructed in %sLoc_db, so linked to things of that one. 
		$ah->{'rawN2rawIdx'}{$sLoc_rawLocS} = $ar_rawIdx_INsLocDb; 
		# The ->{'sorted_rawIdx'} has not been concerned yet. 
	}

	return; 
}# _fill_struct_SEloc () 

=head1 _fill_sorted_rawIdx( \%db_num )

Function : 

           Add ->{'sorted_rawIdx'} => [ [$rawIdx3, $rawIdx8, ...], [$rawIdx2], ... ];
           # Here $rawIdx3_rawNumber is equal to $rawIdx8_rawNumber, and smaller than $rawIdx2_rawNumber;
           # The array_reference in values is linked to ->{'rawN2rawIdx'}{$rawN} for easy changes;

           Add ->{'rawN2srtIdx'}{$rawN} => $sorted_idx_in_'sorted_rawIdx_array';

=cut
sub _fill_sorted_rawIdx {
	my ($ah) = @_; 

	my $tmp_ar = $ah->{'struct'}; 
	for my $ar ( @{ &mathSunhh::_lsAll_inStruct( $ah->{'struct'}, $ah->{'maxP'} ) } ) {
		my ( $t_path, $t_value ) = @$ar; 
		push( @{$ah->{'sorted_rawIdx'}}, $t_value ); # Link this to ->{'rawN2rawIdx'}{$rawN}; 
		$ah->{'rawN2srtIdx'}{ $ah->{'rawNum'}[ $t_value->[0] ] } = $#{ $ah->{'sorted_rawIdx'} }; 
	}
	return; 
} # _fill_srtRawIdx() 

=head1 _hasPos_inLocDb ( \%db_loc, $position )

Return       : (\@realIdx_with_inputPos) === ( [ $realIdx1, $realIdx2, ...  ] )

Description  : This can be further used to extract realLoci by @{$ah->{'realLoc'}}[ @realIdx_with_inputPos ] ; 

=cut
sub _hasPos_inLocDb {
	my ($ah, $posi) = @_; 
	my $rawNum_toChk_ar = mathSunhh::map_windows( 'position' => $posi, 'wind_hash' => $ah->{'wind'} ); 
	my @back; 
	my %has; 
	for my $rawNum_toChk ( @$rawNum_toChk_ar ) {
		defined $ah->{'rawN2rawIdx'}{$rawNum_toChk} or next; 
		# my $realIdx_ar = &mathSunhh::_locDb_rawIdx2realIdx( $ah, $ah->{'rawN2rawIdx'}{$rawNum_toChk}, 1 ); 
		for my $t ( @{ &mathSunhh::_locDb_rawIdx2realIdx( $ah, $ah->{'rawN2rawIdx'}{$rawNum_toChk}, 1 ) } ) {
			defined $has{ $t } and next; 
			$has{ $t } = 1; 
			$posi >= $ah->{'realLoc'}[$t][0] and $posi <= $ah->{'realLoc'}[$t][1] and push(@back, $t); 
		}
	}

	return(\@back); 
}# _hasPos_inLocDb ()


=head1 _map_loc_to_rawIdx ( \%db_loc, [$chkLoc_S, $chkLoc_E] )

Return        : ([$idx1_in_'rawLoc_array', $idx2_in_'rawLoc_array', ...]) 

Description   : 

 Better ignore this inner function. 
 $idx1_in_'rawLoc_array' can be used in $db_loc{'rawLoc'}[$idx1_in_rawLoc_array] (Value=[$rawLocS, $rawLocE]). 
 Be aware that when there are overlaps between rawLocs, the return might be not complete!!! 

=cut
sub _map_loc_to_rawIdx {
	my ($locDB_hr, $locSE_ar) = @_; 
	my @locSE = (@$locSE_ar); 
	$locSE[0] > $locSE[1] and @locSE = @locSE[1,0]; 
	my @back; 

	my $find_idxS = &mathSunhh::_map_number_to_srtIdx( $locDB_hr, $locSE[0], 'both' ); 
	# Check reverse 
	CHK_REV: 
	for ( my $srtIdx=$find_idxS; $srtIdx>=0; $srtIdx-- ) {
		for my $rawIdx (reverse @{ $locDB_hr->{'sorted_rawIdx'}[$srtIdx] }) {
			my ($curLocS, $curLocE) = @{ $locDB_hr->{'rawLoc'}[$rawIdx] }; 
			$curLocE < $locSE[0] and last CHK_REV; 
			$curLocS > $locSE[1] and next; 
			push(@back, $rawIdx); 
		}
	}
	@back = reverse(@back); 
	# Check forward 
	CHK_FWD: 
	for ( my $srtIdx=$find_idxS+1; $srtIdx < @{$locDB_hr->{'sorted_rawIdx'}}; $srtIdx++ ) {
		for my $rawIdx ( @{ $locDB_hr->{'sorted_rawIdx'}[$srtIdx] } ) {
			my ($curLocS, $curLocE)  = @{ $locDB_hr->{'rawLoc'}[$rawIdx] }; 
			$curLocS > $locSE[1] and last CHK_FWD; 
			$curLocE < $locSE[0] and next; 
			push(@back, $rawIdx); 
		}
	}
	return(\@back); 
}# _map_loc_to_rawIdx () 

=head1 _binRealLoc2BinLoc( \%db_loc )

Return       : () 

Function     : 

    When 'binSize' <= 0 , no BIN action is processed. In this case, 'rawLoc' is identical to 'realLoc'. 

    Add ->{'rawIdx2realIdx'}{$rawIdx}    -> [$realIdx_1, $realIdx_2, ...]
          # This relationship is the most important to transform 'rawLoc' to 'realLoc' information. 
    Add ->{'rawLoc'}              => [ [$rawLoc_1_S, $rawLoc_1_E], [$rawLoc_2_S, $rawLoc_2_E], ... ]
          # Here 'rawLoc' stands for 'binLoc'; 
    Add ->{'rawNum'}              => [ $rawLoc_1_S, $rawLoc_2_S, $rawLoc_3_S, ... ]
          # Here 'rawNum' is related to 'rawLoc' instead of 'realLoc'
    Add ->{'wind'}                => mathSunhh::setup_windows(), in which ->{'wind'}{'loci'}
=cut
sub _binRealLoc2BinLoc {
	my ($ah) = @_; 

	if ( $ah->{'binSize'} <= 0 ) {
		@{$ah->{'rawLoc'}} = @{$ah->{'realLoc'}}; 
		@{$ah->{'rawNum'}} = map { $_->[0] } @{$ah->{'rawLoc'}}; 
		for (my $i=0; $i<@{$ah->{'realLoc'}}; $i++) {
			$ah->{'rawIdx2realIdx'}{$i} = [$i]; 
		}
		return; 
	}

	# Sort realNumbers 
	my @realNum = @{ $ah->{'realLoc'} }; 
	my %locDb = %{ &mathSunhh::index_SEloc( \@realNum, { 'binSize' => 0 } ) }; 

	# Construct windows (bins)
	my %wind = %{ mathSunhh::setup_windows( 'ttl_start'=>$ah->{'minV'}, 'ttl_end'=>$ah->{'maxV'}, 'wind_size'=>$ah->{'binSize'}, 'wind_step'=>$ah->{'binSize'}, 'minRatio'=>0 ) }; 
	my %wind_si2wi; 
	for (my $i=0; $i<@{$wind{'info'}{'windSloci'}}; $i++) {
		$wind_si2wi{ $wind{'info'}{'windSloci'}[$i] } = $i; 
	}
	my %used_wi; 
	my %wi2realIdx; # {$wi} = [$realIdx_1, $realIdx_2, ...] 
	for my $tmp_idx_ar (@{$locDb{'sorted_rawIdx'}}) {
		for my $realIdx ( @$tmp_idx_ar ) {
			# In this way, all $realIdx are sorted according to their original SEloc. 
			# So, the order in 'rawIdx2realIdx' is also sorted according to their original SEloc.
			my @real_SE = @{$ah->{'realLoc'}[$realIdx]}; 
			my @idxS_arr = @{ mathSunhh::map_windows( 'position' => $real_SE[0], 'wind_hash'=>\%wind ) }; 
			my @idxE_arr = @{ mathSunhh::map_windows( 'position' => $real_SE[1], 'wind_hash'=>\%wind ) }; 
			my ($si_S, $si_E) = &mathSunhh::minmax( [@idxS_arr, @idxE_arr] ); 
			my ($wi_S, $wi_E) = ( $wind_si2wi{$si_S}, $wind_si2wi{$si_E} ); 
			for my $t_wi ($wi_S .. $wi_E) {
				$used_wi{$t_wi} = 1; 
				push(@{ $wi2realIdx{$t_wi} }, $realIdx); 
			}
		}
	}
	my @used_wi_arr = sort {$a<=>$b} keys %used_wi; 
	for my $t_wi (@used_wi_arr) {
		my $t_si = $wind{'info'}{'windSloci'}[$t_wi]; 
		push(@{$ah->{'rawLoc'}}, [ @{$wind{'loci'}{$t_si}}[0,1] ]); 
		$ah->{'rawIdx2realIdx'}{ $#{$ah->{'rawLoc'}} } = [ @{ $wi2realIdx{$t_wi} } ]; 
		push(@{$ah->{'rawNum'}}, $ah->{'rawLoc'}[-1][0]); 
	}
	$ah->{'wind'} = \%wind; 
	
	return; 
}# _binRealLoc2BinLoc () 

=head1 _closest_number( $number_toCheck, \@source_numbers )

Return       : ( $number_from_source_closest_to_number_toCheck )

Description  : If there are two numbers with same distance flanking $number_toCheck, the smaller one will be chosen. 

=cut
sub _closest_number {
	my ($chkN, $ar) = @_; 
	@$ar > 0 or &stopErr("[Err] No database numbers for checking in\n"); 
	my %best = ( 'idx' => 0, 'abs_diff' => abs( $ar->[0]-$chkN ) ); 
	my $cnt = -1; 
	for my $dbN (@$ar) {
		$cnt ++; 
		my $diff = abs( $dbN - $chkN ); 
		if ($diff < $best{'abs_diff'}) {
			$best{'idx'} = $cnt; 
			$best{'abs_diff'} = $diff; 
		}
	}

	return( $ar->[$best{'idx'}] ); 
}# _closest_number ()

1; 
