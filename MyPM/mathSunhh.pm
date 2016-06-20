package mathSunhh; 
#BEGIN {
#	push(@INC,'/usr/local/share/perl5/'); 
#}
# This is a package storing sub functions. 
# So there isn't objects created. 
use strict; 
use warnings; 
use Statistics::Descriptive; 
use Scalar::Util qw(looks_like_number blessed);
use Exporter qw(import);
our @EXPORT = qw(ins_calc);
our @EXPORT_OK = qw();

use LogInforSunhh; 

my $ms_obj = mathSunhh->new(); 

############################################################
#  Methods
############################################################

sub new {
	my $class = shift; 
	my $self = {}; 
	bless $self, $class; 
	
	$self->_initialize(@_); 
	
	return $self; 
}

sub _initialize {
	my $self = shift; 
	my %parm = $self->_setHashFromArr(@_); 
	return; 
}


=head2 repArr (\@units, 'each'=>1, 'times'=>1, 'length'=>\@units) 

Required : \@units

Function : A repeat function similar to rep() in R. 
           First apply 'each', then apply 'times', and at least apply 'length'; 

Return   : \@repeated_eles 

=cut
sub repArr {
	my $self = shift; 
	my $ref_arr = shift; 
	my %parm = $self->_setHashFromArr(@_); 
	$parm{'each'} = $parm{'each'} // 1; 
	$parm{'times'} = $parm{'times'} // 1; 
	$parm{'length'} = $parm{'length'} // ( ($#$ref_arr+1) * $parm{'each'} * $parm{'times'} ); 
	$parm{'length'} > 0 or return []; 
	my @each_arr; 
	for my $te (@$ref_arr) {
		push( @each_arr, ( $te ) x $parm{'each'} ); 
	}
	my @time_arr; 
	for ( my $i=0; $i<$parm{'times'}; $i++ ) {
		push(@time_arr, @each_arr); 
	}
	my @back_arr; 
	for ( my ($i, $time_idx) = (0,0) ; $i<$parm{'length'}; $i++ ) {
		$time_idx > $#time_arr and $time_idx -= ( $#time_arr + 1 ); 
		push(@back_arr, $time_arr[$time_idx]); 
		$time_idx ++; 
	}
	return \@back_arr; 
}# End repArr() 


=head2 newNumber ( 'other_safeNumber'=>[$mathSunhhObj->{'safeNumber'}, ...], 'onlyMerge'=>0, 'debug'=>0 )

Required : Null 

Function : Construct 'safeNumber' hash to store not used numbers in the same object. 
            It seems that this number is not only unique in the same object, but also unique in all objects in the same module calling this mathSunhh.pm module. 

Input    : Null 

Return   : A new number not used before. 

=cut
sub newNumber {
	my $self = shift; 
	my %parm = $self->_setHashFromArr(@_); 
	$parm{'debug'} = $parm{'debug'} // 0; 
	$parm{'onlyMerge'} = $parm{'onlyMerge'} // 0; 
	if ( defined $parm{'other_safeNumber'} ) {
		for my $h1 ( @{$parm{'other_safeNumber'}} ) {
			my ($min1, $max1) ; 
			for my $tN ( keys %{ $h1->{'has'} } ) {
				defined $min1 or $min1 = $tN; 
				defined $max1 or $max1 = $tN; 
				$min1 > $tN and $min1 = $tN; 
				$max1 < $tN and $max1 = $tN; 
				$parm{'debug'} and defined $self->{'safeNumber'} and defined $self->{'safeNumber'}{'has'}{$tN} and &tsmsg("[Wrn] existed number $tN in other_safeNumber\n"); 
				$self->{'safeNumber'}{'has'}{$tN} = 1; 
			}
			defined $self->{'safeNumber'}{'min'} or $self->{'safeNumber'}{'min'} = $min1; 
			defined $self->{'safeNumber'}{'max'} or $self->{'safeNumber'}{'max'} = $max1; 
			defined $self->{'safeNumber'}{'has'} or $self->{'safeNumber'}{'has'} = {}; 
			defined $min1 and $min1 < $self->{'safeNumber'}{'min'} and $self->{'safeNumber'}{'min'} = $min1; 
			defined $max1 and $max1 > $self->{'safeNumber'}{'max'} and $self->{'safeNumber'}{'max'} = $max1; 
		}
	}
	$parm{'onlyMerge'} and return ; 
	if ( defined $self->{'safeNumber'} ) {
		my $newN = ( defined $self->{'safeNumber'}{'max'} ) ? $self->{'safeNumber'}{'max'} + 1 : 1; 
		while ( defined $self->{'safeNumber'}{'has'}{$newN} ) {
			$newN ++; 
		}
		$self->{'safeNumber'}{'max'} = $newN; 
		$self->{'safeNumber'}{'has'}{$newN} = 1; 
		$self->{'safeNumber'}{'lastN'} = $newN; 
		return $newN; 
	} else {
		my $newN = 1; 
		$self->{'safeNumber'}{'min'} = $newN; 
		$self->{'safeNumber'}{'max'} = $newN; 
		$self->{'safeNumber'}{'has'}{$newN} = 1; 
		$self->{'safeNumber'}{'lastN'} = $newN; 
		return $newN; 
	}
}# newNumber() 

=head2 offspringArray( $rootID, $code_reference_of_function, 'unique'=>1 )

Required : 
  $rootID
  $code_reference_of_function : sub { return [@arr_value_as_next_self_loop] }

Function : Trace back all offsprings from root ID according to sub_routine_reference given. 
           Be aware that if the offspring has the same ID of rootID, this function will terminate!!! 
           This is used to avoid infinite cycles caused by relationship: rootID is a child of rootID. 

Input    : ( $rootID, sub { return @offspring; }, 'unique'=>1 )

Output   : [offspring1, offspring2, ...]

=cut
sub offspringArray {
	my $self = shift;
	my $rootID = shift;  # rootID from which to find offsprings. 
	my $coderef = shift; # reference of subroutine which return a array of offspring list. 
	my %parm = $self->_setHashFromArr(@_); 
	$parm{'unique'} = $parm{'unique'} // 1; 

	my @cIDs;            # All Children IDs. 
	my @off1 = ( &$coderef( $rootID ) ); # Offsprings directly from rootID. 
	scalar(@off1) == 0 and return [];
	my %haveID ; 
	$haveID{$rootID} = 1; 

	for my $cID ( @off1 ) {
		defined $haveID{$cID} and next; 
		push(@cIDs, $cID);
		my @ccIDs = grep { !(defined $haveID{$_}) } @{ $self->offspringArray( $cID, $coderef, %parm ) };
		push(@cIDs, @ccIDs);
	}

	if ( $parm{'unique'} ) {
		my %usedID; 
		my @bIDs; 
		for my $cID (@cIDs) {
			defined $usedID{$cID} or push(@bIDs, $cID); 
			$usedID{$cID} = 1; 
		}
		@cIDs = @bIDs; 
	}

	return [@cIDs];
}# _allChild()

=head2 setup_windows( 'ttl_start'=>1, 'ttl_end'=>99999999, 'wind_size'=>1000, 'wind_step'=> // 'wind_size', 
  'minRatio'=> 0
)

Required: 
 

Function: return a list of windows in hash reference according to given [window_size, window_step, total_start, total_end]

Return  : \%back_wind; 
 'keys'  => values
 'info' => {'ttl_start/ttl_end/wind_size/wind_step/minRatio/windSloci' => values} 
 'loci'  => {
              'wind_start_position' => [start_pos, end_pos, interval_len]
            }
 
=cut
sub setup_windows {
	my $self = shift; 
	my %parm = $self->_setHashFromArr(@_); 
	$parm{'ttl_start'} = $parm{'ttl_start'} // 1; 
	$parm{'ttl_end'}   = $parm{'ttl_end'}   // 99999999; 
	$parm{'wind_size'} = $parm{'wind_size'} // 1000; 
	$parm{'wind_step'} = $parm{'wind_step'} // $parm{'wind_size'}; 
	$parm{'minRatio'}  = $parm{'minRatio'}  // 0; 
	my %back_wind; 
	for (qw/ttl_start ttl_end wind_size wind_step minRatio/) {
		$back_wind{'info'}{$_} = $parm{$_}; 
	}
	my $min_windSize = int($parm{'minRatio'}*$parm{'wind_size'}); 
	$min_windSize < $parm{'minRatio'}*$parm{'wind_size'} and $min_windSize ++; 
	$min_windSize == 0 and $min_windSize = 1;  
	
	for (my $si=$parm{'ttl_start'}; $si+$min_windSize-1 <= $parm{'ttl_end'}; $si += $parm{'wind_step'}) {
		my $ei = $si + $parm{'wind_size'}-1; 
		$ei > $parm{'ttl_end'} and $ei = $parm{'ttl_end'}; 
		my $cur_len = $ei-$si+1; 
		$back_wind{'loci'}{$si} = [$si, $ei, $cur_len]; 
		push(@{$back_wind{'info'}{'windSloci'}}, $si); 
	}
	return \%back_wind; 
}# sub setup_windows

=head2 map_windows ( 'posi/position'=>Integer, 'wind_hash'=>setup_windows->(), 
  'ttl_start'=>1, 'ttl_end'=>99999999, 'wind_size'=>1000, 'wind_step'=>'wind_size', 'minRatio'=>0
)

Required: 
 'posi/position'

Function: 
 Given a position, return an array reference recording all start_positions of windows that this position locates in. 
 'wind_hash' will mask 'ttl_*' and 'wind_size/wind_step'. 

Return  : 
 \@back_si = [si_1, si_2, si_3, ...]
 Here "si" should be a key of %{$parm{'wind_hash'}{loci}}; 

=cut
sub map_windows {
	my $self = shift; 
	my %parm = $self->_setHashFromArr(@_); 
	my $posi = $parm{'posi'} // $parm{'position'} // &stopErr("[Err] No position assigned.\n"); 
	Scalar::Util::looks_like_number( $posi ) or &stopErr("[Err] input position [$posi] is not like a number.\n"); 
	$posi = int($posi); 
	if (defined $parm{'wind_hash'}) {
		for ( qw/ttl_start ttl_end wind_size wind_step minRatio/ ) {
			$parm{$_} = $parm{$_} // $parm{'wind_hash'}{'info'}{$_} // undef(); 
		}
	}
	$parm{'ttl_start'} = $parm{'ttl_start'} // 1; 
	$parm{'ttl_end'}   = $parm{'ttl_end'}   // 99999999; 
	$parm{'wind_size'} = $parm{'wind_size'} // 1000; 
	$parm{'wind_step'} = $parm{'wind_step'} // $parm{'wind_size'}; 
	$parm{'minRatio'}  = $parm{'minRatio'}  // 0; 
	
	my @back_si; 
	my $resid = ($posi-$parm{'ttl_start'}) % $parm{'wind_step'}; 
	my $end_si = $posi - $resid; 
	
	my $min_windSize = int($parm{'minRatio'}*$parm{'wind_size'}); 
	$min_windSize < $parm{'minRatio'}*$parm{'wind_size'} and $min_windSize ++; 
	$min_windSize == 0 and $min_windSize = 1;  
	
	for (my $si=$end_si; $si+$parm{'wind_size'}-1>=$posi && $si >= $parm{'ttl_start'}; $si-=$parm{'wind_step'}) {
		$si+$min_windSize-1 <= $parm{'ttl_end'} or next; 
		push(@back_si, $si); 
	}

	@back_si = reverse(@back_si); 
	
	return \@back_si; 
}# sub map_windows 


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
#		return &min($e1, $e2) - &max($s1, $s2) + 1; 
#	}
#}
=head2 ovl_region( $start1, $end1, $start2, $end2 )

Required: ($start1, $end1, $start2, $end2)

Function: Check overlapping region

Return  : ( $overlapped_bp_length, [ ovl_start, ovl_end ] )

=cut
sub ovl_region {
	my $self = shift; 
	my ($s1, $e1, $s2, $e2) = @_; 
	$s1 > $e1 and ($s1, $e1) = ($e1, $s1); 
	$s2 > $e2 and ($s2, $e2) = ($e2, $s2); 
	if ( $e1 < $s2 or $s1 > $e2 ) {
		return (0, []); 
	} else {
		my $ee = &min($e1, $e2); 
		my $ss = &max($s1, $s2); 
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
	my $self = shift; 
	my $list1 = shift; 
	my $list2 = shift; 
	my %parm = $self->_setHashFromArr(@_); 
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
				my ($ovl_len, $ovl_loc) = $self->ovl_region($s1, $e1, $s2, $e2); 
				$ovl_len > 0 and do { $ovlCnt ++; push( @ovlLoc, [@$ovl_loc] );  }; 
			}
		}
		# Merge @ovlLoc blocks
		@ovlLoc = @{ $self->mergeLocBlk(\@ovlLoc) }; 
		
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
					my ( $ovl_o1, $ovl_o2 ) = $self->ovl_region($s1, $e1, $s3, $e3); 
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
					my ( $ovl_o1, $ovl_o2 ) = $self->ovl_region($s1, $e1, $s3, $e3); 
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

 Can be used as &mergeLocBlk( [ [$s11,$e11], [$s12,$e12] ], [[$s21, $e21]], 'dist2join'=>1  )
 
 Here 'dist2join' is the number distance between two blocks that can be joined. Default is 0, means [4,6] and [7,10] will be two blocks. 
 If set 'dist2join'=>2, [4,6] and [8,10] will be joined into [4,10] ; 

Return     : [[ss1, ee1], [ss2, ee2], ...]

=cut 
sub mergeLocBlk {
	my $self = shift; 
	my @back_blk; 
	my @srt_blk; 
	my @para_array; 
	for my $blk_arr (@_) {
		ref($blk_arr) eq '' and do { push(@para_array, $blk_arr); next; }; 
		push(@srt_blk, @$blk_arr); 
	}
	my %parm = &_setHashFromArr(@para_array); 
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

=head1 _setHashFromArr(@keyVal_array)

=head2 _setHashFromArr(@keyVal_array)

Required: @keyVal_array

Function: @keyVal_array contain ( key1, val1, key2, val2, ... ) pairs. 
          It will skip the latter duplicated keyN. 

Return  : %back_hash
  In %back_hash : 
   {keyN} => valN

=cut
sub _setHashFromArr {
	my $self = shift; 
	if ( ref($self) ne 'mathSunhh' ) {
		ref($self) eq '' or &stopErr("[Err] _setHashFromArr() input [$self] illegal.\n"); 
		defined $self and unshift(@_, $self); 
	}
	my %back_hash; 
	for (my $i=0; $i<@_; $i+=2) {
		my $val; 
		if (exists $_[$i+1]) {
			$val = $_[$i+1]; 
		} else {
			exists $back_hash{$_[$i]} or &tsmsg("[Wrn] Input array is not even! Use undef() for key [", $_[$i],"]\n"); 
			$val = undef(); 
		}
		exists $back_hash{$_[$i]} or $back_hash{$_[$i]} = $val; 
	}
	return(%back_hash); 
}# _setHashFromArr() 

=head2 switch_position('qry2ref'=>\%refPos, 'qryID'=>$queryID, 'qryPos'=>$queryPos, 'qryStr'=>'+/-')

Input       : 
 %refPos   : {qryID} => sortedByQrySE [ [qryS_1, qryE_1, refID_1, refS_1, refE_1, refStr_1(+/-)], [qryS_2, qryE_2, refID_2, refS_2, refE_2, refStr_2], ... ]
 
               qryS <= qryE && refS <= refE 

			   qryS_1 <= qryS_2 ... 

Function    : 

 Get position on ref for qryID_qryPos

Return      : 

 ([refID_1, refPos_1, refStr_1], [refID_2, refPos_2, refStr_2])

=cut
sub switch_position {
	my $self = shift; 
	my %parm = $self->_setHashFromArr(@_); 
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

############################################################
#  Sub-routines. 
############################################################


=head1 get_xy_byScale( 'start_xy'=>[$sX,$sY], 'end_xy'=>[$eX,$eY], 'se_position'=>[$sPos,$ePos], 'need_position'=>[$p1,$p2,...] )

Return : ( [$needP_1_x, $needP_1_y], [$needP_2_x, $needP_2_y], ... )

=cut
sub get_xy_byScale {
	my %parm = &_setHashFromArr(@_); 
	my @back; 
	
	$parm{'se_position'}[1] == $parm{'se_position'}[0] and do { @back = ([@{$parm{'start_xy'}}]); return(@back); }; 
	
	my $scale_x = ($parm{'end_xy'}[0]-$parm{'start_xy'}[0])/($parm{'se_position'}[1]-$parm{'se_position'}[0]); 
	my $scale_y = ($parm{'end_xy'}[1]-$parm{'start_xy'}[1])/($parm{'se_position'}[1]-$parm{'se_position'}[0]); 
	for my $curP (@{$parm{'need_position'}}) {
		my $curP_x = $parm{'start_xy'}[0] + $scale_x * ($curP-$parm{'se_position'}[0]) ; 
		my $curP_y = $parm{'start_xy'}[1] + $scale_y * ($curP-$parm{'se_position'}[0]) ; 
		push(@back, [$curP_x, $curP_y]); 
	}
	
	return(@back); 
}# get_xy_byScale()

=head1 cnvt_to_rgb ( $min, $max, $val, $col ) 

Function : Get a series of color for heat map. 

Return : ( "rgb(R,G,B)" )

Required : 
  $col = [ [r1,g1,b1], [r2,g2,b2], [r3,g3,b3], ... ]
  $min , $max : value range; 
  $val : current value. 

Reference: http://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map

=cut
sub cnvt_to_rgb {
	my ( $min, $max, $val, $col ) = @_; 
	my $col_i = $#$col;
	my $v = ($val-$min)/($max-$min) * $col_i;
	my $i1 = &min( int($v), $col_i );
	my $i2 = &min( int($v)+1, $col_i);
	my $f = $v - $i1;
	my @back;
	for (my $j=0; $j<3; $j++) {
		$back[$j] = int( $col->[$i1][$j] + $f * ($col->[$i2][$j]-$col->[$i1][$j]) );
	}
	return( sprintf("rgb(%.0f,%.0f,%.0f)", $back[0], $back[1], $back[2])  );
}# cnvt_to_rgb () 



=head1 _parseCol( $col_string )

Required    : $col_string
 $col_string : Like 0-100 | 0,3,4,5-7,1-10 | 10-8

Function    : 
 '0-5' returns (0 .. 5)
 '0,3,5-7,6-9' returns (0,3,5,6,7,6,7,8,9)
 ' -3 - -4, 10-8'  returns (-3,-4,10,9,8)

Return      : (@col_numbers)

=cut
sub _parseCol {
	my @ncols; 
	for my $tc (split(/,/, shift)) {
		$tc =~ s!\s!!g; 
		$tc =~ m!^$! and next; 
		if ( $tc =~ m!^\-?\d+$! ) {
			push(@ncols, $tc); 
		} elsif ( $tc =~ m!^(\-?\d+)\-(\-?\d+)$! ) {
			my ($s, $e) = ($1, $2); 
			if ( $s <= $e ) {
				push(@ncols, ($s .. $e)); 
			} else {
				push(@ncols, reverse($e .. $s)); 
			}
		} else {
			&stopErr("[Err] Unparsable column tag [$tc] in _parseCol().\n"); 
		}
	}
	return (@ncols); 
}# _parseCol() 

=head1 _mean(@numbers)
=cut
sub _mean {
	my $stat = Statistics::Descriptive::Full->new(); 
	$stat->add_data(@_); 
	return $stat->mean(); 
}
=head1 _sum(@numbers)
=cut
sub _sum {
	my $stat = Statistics::Descriptive::Full->new(); 
	$stat->add_data(@_); 
	return $stat->sum(); 
}
=head1 _median(@numbers)
=cut
sub _median {
	my $stat = Statistics::Descriptive::Full->new(); 
	$stat->add_data(@_); 
	return $stat->median(); 
}


=head1 min(@numbers)

=head2 min(@numbers)

Function: This is not a method, but a sub-routine()

=cut
sub min {
	my $min = shift; 
	unless ( ref($min) eq '' ) {
		ref($min) eq 'mathSunhh' or &stopErr("[Err] min() input should be an array of number.\n"); 
		$min = shift; 
	}
	for (@_) {
		defined $_ or next; 
		defined $min or $min = $_; 
		$min > $_ and $min = $_; 
	}
	return $min; 
}# min() 

=head1 max(@numbers)

=head2 max(@numbers)

Function: This is not a method, but a sub-routine()

=cut
sub max {
	my $max = shift; 
	unless ( ref($max) eq '' ) {
		ref($max) eq 'mathSunhh' or &stopErr("[Err] max() input should be an array of number.\n"); 
		$max = shift; 
	}
	for (@_) {
		defined $_ or next; 
		defined $max or $max = $_; 
		$max < $_ and $max = $_; 
	}
	return $max; 
}# max() 

=head1 ins_calc( \@numbers, $min_valid_number_count )

Function: This is not a method, but a sub-routine(). 

Description: For calculating insert sizes. 
             Following Heng Li's bwa method (Estimating Insert Size Distribution). 
             But the max/min distance of INS are only using 6 * sigma values. 
             http://linux.die.net/man/1/bwa
             BWA estimates the insert size distribution per 256*1024 read pairs. It first collects pairs of reads with both ends mapped with a single-end quality 20 or higher and then calculates median (Q2), lower and higher quartile (Q1 and Q3). It estimates the mean and the variance of the insert size distribution from pairs whose insert sizes are within interval [Q1-2(Q3-Q1), Q3+2(Q3-Q1)]. The maximum distance x for a pair considered to be properly paired (SAM flag 0x2) is calculated by solving equation Phi((x-mu)/sigma)=x/L*p0, where mu is the mean, sigma is the standard error of the insert size distribution, L is the length of the genome, p0 is prior of anomalous pair and Phi() is the standard cumulative distribution function. For mapping Illumina short-insert reads to the human genome, x is about 6-7 sigma away from the mean. Quartiles, mean, variance and x will be printed to the standard error output.

Input      : (\@ins_value_array)

Output     : (\%hash_of_values) 
             keys = qw(SUM COUNT MEAN MEDIAN Q1 Q3 interval_low interval_high interval_mean interval_median interval_var interval_stdev interval_cnt limit_low limit_high min max)

=cut
sub ins_calc {
	my $r_arr = shift; 
	unless ( ref($r_arr) eq 'ARRAY' ) {
		ref($r_arr) eq 'mathSunhh' or &stopErr("[Err] ins_calc() input-1st should be an array reference.\n"); 
		$r_arr = shift; 
	}
	my $min_val_number = shift // 1; 
	my %back; 
	if ( (! defined $r_arr) or scalar(@$r_arr) < $min_val_number) {
		for my $ta (qw/SUM COUNT MEAN MEDIAN Q1 Q3 interval_low interval_high interval_mean interval_median interval_var interval_stdev limit_low limit_high/) {
			$back{$ta} = ''; 
		}
		return \%back; 
	}
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@$r_arr); 
	$back{'SUM'} = $stat->sum(); 
	$back{'COUNT'} = $stat->count(); 
	$back{'MEAN'} = $stat->mean(); 
	$back{'MEDIAN'} = $stat->median(); 
	$back{'Q1'} = $stat->quantile(1); 
	$back{'Q3'} = $stat->quantile(3); 
	$back{'interval_low'}  = $back{'Q1'} - 2 * ($back{'Q3'}-$back{'Q1'}); 
	$back{'interval_high'} = $back{'Q3'} + 2 * ($back{'Q3'}-$back{'Q1'}); 
	$back{'min'} = $stat->min(); 
	$back{'max'} = $stat->max(); 
	
	$stat->clear(); 
	my @sub_arr; 
	my $ti = -1; 
	for my $ta (@$r_arr) {
		$ti ++; 
		$ta >= $back{'interval_low'} and $ta <= $back{'interval_high'} and do { push(@sub_arr, $ta); push(@{$back{'interval_idx'}}, $ti); }; 
	}
	$stat->add_data(@sub_arr); 
	$back{'interval_cnt'}   = scalar( @sub_arr ); 
	$back{'interval_arr'}   = \@sub_arr; 
	$back{'interval_mean'}  = $stat->mean(); 
	$back{'interval_median'} = $stat->median(); 
	$back{'interval_var'}   = $stat->variance(); 
	$back{'interval_stdev'} = $stat->standard_deviation(); 
	$back{'limit_low'}  = $back{'interval_mean'} - 6 * $back{'interval_stdev'}; 
	$back{'limit_high'} = $back{'interval_mean'} + 6 * $back{'interval_stdev'}; 
	$stat->clear(); 
	return \%back; 
}# ins_calc() 

=head1 permutations( \@list_of_element, $number_of_list ) 

This is not a method, but a sub-routine. 

Given (\@list_of_ele, $n_in_class), return all permutations by array. Return ([@perm1_of_ele], [@perm2], ...)

=cut
sub permutations {
	my $list = shift; 
	unless ( ref($list) eq 'ARRAY' ) {
		ref($list) eq 'mathSunhh' or &stopErr("[Err] permutations() input-1st should be an array reference.\n"); 
		$list = shift; 
	}
	my $n = shift; 
	$n = $n // scalar(@$list); 
	$n > @$list and return ($list); 
	$n <= 1 and return(map {[$_]} @$list); 
	my @perm; 
	for my $i (0 .. $#$list) {
		my @rest = @$list; 
		my $val = splice(@rest, $i, 1); 
		for ( &permutations(\@rest, $n-1) ) {
			push(@perm, [$val, @$_]); 
		}
	}
	return @perm; 
}#sub permutations() 

=head1 combinations(\@list_of_ele, $n_in_class)

This is not a method, but a sub-routine. 

Given (\@list_of_ele, $n_in_class), return all combinations by array. Return ([@comb1_of_ele], [@comb2_of_ele], ...) 

=cut
sub combinations {
	my $list = shift; 
	unless ( ref($list) eq 'ARRAY' ) {
		ref($list) eq 'mathSunhh' or &stopErr("[Err] combinations() input-1st should be an array reference.\n"); 
		$list = shift; 
	}
	my $n = shift; 
	$n = $n // scalar(@$list); 
	$n > @$list and return ($list); 
	$n <= 1 and return(map {[$_]} @$list); 
	my @comb; 
	for (my $i=0; $i+$n<=@$list; $i++) {
		my $val = $list->[$i]; 
		my @rest = @$list[$i+1 .. $#$list]; 
		for (&combinations(\@rest, $n-1)) {
			push(@comb, [$val, @$_]); 
		}
	}
	return @comb; 
}#sub combinations

=head1 dvd_array( \@array_to_be_divided, $number_of_subgroups, $If_keep_order[default=0], $prev_grp_colN )

Function : Divide @array_to_be_divided into $number_of_subgroups subgroups. 
If $If_keep_order is TRUE, the elements in subgroups will be sequenctial same to the raw order. 

Return   : ( [ \@subgroup_1, \@subgroup_2, ... ] )

Example  : 
  &dvd_array( [0,1,2,3,4,5,6,7] , 3 )    returns [ [0,3,6], [1,4,7], [2,5] ] ; 
  &dvd_array( [0,1,2,3,4,5,6,7] , 3, 1 ) returns [ [0,1,2], [3,4,5], [6,7] ] ; 

=cut
sub dvd_array {
	my $inA = shift; 
	my $grpN = shift; 
	my $keepOrder = shift; 
	my $prevGrp_ColN = shift; 
	$prevGrp_ColN //= 'N'; 
	$grpN = int($grpN);
	$grpN >= 1 or $grpN = 1; 
	my @sub_grps;
	if ( $keepOrder ) {
		my $num_ttl = scalar(@$inA); 
		my $num_subG = int( $num_ttl/$grpN ); 
		$num_subG * $grpN >= $num_ttl or $num_subG++; 
		$num_subG * $grpN >= $num_ttl or &stopErr("[Err] num_ttl=$num_ttl , grpN=$grpN\n"); 
		if ( $prevGrp_ColN eq 'N' ) {
			for (my $i=0; $i<@$inA; $i+=$num_subG) {
				my $e = $i+$num_subG-1; 
				$e > $#$inA and $e = $#$inA; 
				push( @sub_grps, [ @{$inA}[ $i .. $e ] ] ); 
			}
		} else {
			my $prev_tag; 
			my $curr_num = 0; 
			my $subGrp_idx = 0; 
			my $num_cutoff = $num_subG; 
			for (sort { $a->[$prevGrp_ColN]<=>$b->[$prevGrp_ColN] || $a->[$prevGrp_ColN] cmp $b->[$prevGrp_ColN] } @$inA) {
				$prev_tag //= $_->[$prevGrp_ColN]; 
				if ( $prev_tag ne $_->[$prevGrp_ColN] ) {
					$prev_tag = $_->[$prevGrp_ColN]; 
					if ( $curr_num < $num_cutoff ) {
						push( @{$sub_grps[$subGrp_idx]}, $_ ); 
					} else {
						$subGrp_idx ++; 
						$num_cutoff = $curr_num + $num_subG; 
						push( @{$sub_grps[$subGrp_idx]}, $_ ); 
					}
				} else {
					push( @{$sub_grps[$subGrp_idx]}, $_ ); 
				}
				$curr_num ++; 
			}
		}
	} else {
		$prevGrp_ColN ne 'N' and &stopErr("[Err] mathSunhh::dvd_array(): prev_grp_colN can be only ['N'] when If_keep_order is 0.\n"); 
		for (my $i=0; $i<@$inA; $i+=$grpN) {
			for (my $j=0; $j<$grpN; $j++) {
				my $k = $j+$i;
				$k < @$inA or last;
				push(@{$sub_grps[$j]}, $inA->[$k]);
			}
		}
	}
	return \@sub_grps;
}# dvd_array ()


=head1 complete_pair_for_grouping ( \%relation_pairs )

Input      : Format of %relation_pairs is 
             {$element_1}{$element_2} = 1; 
             {$element_3}{$element_4} = 1; 
             {$element_4}{$element_3} = 1; 
             {$element_5}{$element_6} = 1; 
             {$element_1}{$element_6} = 1; 
             {$element_6}{$element_1} = 1; 
             ...... 

Output     : ( \%complete_relation_pairs )
             {$element_1}{$element_2} = 1; 
             {$element_2}{$element_1} = 1; 
             {$element_3}{$element_4} = 1; 
             {$element_4}{$element_3} = 1; 
             {$element_5}{$element_6} = 1; 
             {$element_5}{$element_5} = 1; 
             {$element_1}{$element_6} = 1; 
             {$element_6}{$element_1} = 1; 
             ...... 

=cut
sub complete_pair_for_grouping {
	my ( $pair_href ) = @_; 
	my %back_hash = %$pair_href; 
	my %to_add; 
	for my $tk1 ( keys %back_hash ) {
		for my $tk2 ( keys %{$back_hash{$tk1}} ) {
			defined $back_hash{$tk2} and defined $back_hash{$tk2}{$tk1} and next; 
			$to_add{$tk2}{$tk1} = 1; 
		}
	}
	for my $tk1 ( keys %to_add ) {
		for my $tk2 ( keys %{$to_add{$tk1}} ) {
			$back_hash{$tk1}{$tk2} = 1; 
		}
	}
	return (\%back_hash); 
} # complete_pair_for_grouping() 

=head1 divide_group_fromHash ( \%relation_pairs ) 

Required   : Do be aware that in %relation_pairs, each pair should have two formed key pairs, {k1}{k2} and {k2}{k1}

Input      : Format of %relation_pairs is 
             {$element_1}{$element_2} = 1; 
             {$element_2}{$element_1} = 1; 
             {$element_3}{$element_4} = 1; 
             {$element_4}{$element_3} = 1; 
             {$element_5}{$element_6} = 1; 
             {$element_5}{$element_5} = 1; 
             {$element_1}{$element_6} = 1; 
             {$element_6}{$element_1} = 1; 
             ...... 

Function   : Divide elements in %relation_pairs into groups. 
             Each group contains related elements, and there is no related pairs between groups. 
             This should be much faster than &divide_group(); 

Output     : (\@back_groups_array)
 @back_groups_array format is ( [group_1_ele_1, group1_ele_2, ...], [group_2_ele_1, group_2_ele_2, ...], ... ); 

=cut
sub divide_group_fromHash {
	my ( $pair_href ) = @_; 
	my @array_input; 
	map { push( @array_input, [ $_, keys %{$pair_href->{$_}} ] ); } keys %$pair_href; 
	my ($back_aref) = &divide_group_fromArray(\@array_input); 
	
	return ($back_aref); 
}# divide_group_fromHash() 

=head1 divide_group_fromArray ( \@elements_in_groups ) 

Required   : 

Input      : Format of @elements_in_groups is 
             ([ ele1, ele2, ... ], 
             [ele3, ele4, ele5, ...], 
             [ele1, ele3], 
             [ele6, ele7, ele8], 
             [ele11, ele10, ele9], 
             ......)

Function   : Divide elements in @elements_in_groups into groups. 
             Each group contains related elements, and there is no related pairs between groups. 

Output     : (\@back_groups_array)
 @back_groups_array format is ( [group_1_ele_1, group1_ele_2, ...], [group_2_ele_1, group_2_ele_2, ...], ... ); 

=cut
sub divide_group_fromArray {
	my ( $pair_aref ) = @_; 
	my (%eleID_to_grpID, %grpID_wi_eleID, $grpID); 
	$grpID = 0; 
	for my $ar1 ( @$pair_aref ) {
		my @ta = @$ar1; 
		my ($minGID, $maxGID) = ('X', 'X'); 
		my %usedGID; 
		# Check if $ele2 in @ta already exists in %eleID_to_grpID, and determine the $minGID and maxGID of parent groups ; 
		for my $ele2 ( @ta ) {
			defined $eleID_to_grpID{$ele2} or next; 
			$usedGID{ $eleID_to_grpID{$ele2} } = 1; 
			$minGID = ( $minGID eq 'X' )
			  ? $eleID_to_grpID{$ele2} 
			  : ( $minGID > $eleID_to_grpID{$ele2} ) 
			    ? $eleID_to_grpID{$ele2} 
			    : $minGID 
			; 
			$maxGID = ( $maxGID eq 'X' ) 
			  ? $eleID_to_grpID{$ele2}
			  : ( $maxGID < $eleID_to_grpID{$ele2}  )
			    ? $eleID_to_grpID{$ele2}
			    : $maxGID
			; 
		}# End for my $ele2 ( @ta ) 
		
		if ( $minGID eq 'X' ) {
			# Setup a new group if there is no existing parent group; 
			map { $eleID_to_grpID{$_} = $grpID; $grpID_wi_eleID{$grpID}{$_}=1; } @ta; 
			$grpID ++; 
		} elsif ( $minGID < $maxGID ) {
			# Move groups with GID larger than $minGID to group with $minGID; 
			# And then remove those larger groups. 
			# Add all of the rest elements (eleID) to group of $minGID; 
			delete $usedGID{ $minGID } ; 
			for my $t_GID ( keys %usedGID ) {
				map { $eleID_to_grpID{$_}=$minGID; $grpID_wi_eleID{$minGID}{$_}=1; } keys %{ $grpID_wi_eleID{$t_GID} }; 
				delete $grpID_wi_eleID{$t_GID}; 
			}
			## Adding rest elements. 
			for my $eleID (@ta) {
				exists $eleID_to_grpID{$eleID} and next; 
				$eleID_to_grpID{$eleID}=$minGID; 
				$grpID_wi_eleID{$minGID}{$eleID}=1; 
			}
		} else {
			# Here minGID == maxGID, which means all of the elements are in the same existing group. 
			# So simply add the ones not included previously to that group. 
			for my $eleID (@ta) {
				exists $eleID_to_grpID{$eleID} and next; 
				$eleID_to_grpID{$eleID}=$minGID; 
				$grpID_wi_eleID{$minGID}{$eleID}=1; 
			}
		}#End if ( $minGID eq 'X' ) 
	}#End for my $ar1 ( @$pair_aref ) 
	
	my @back_array; 
	my (%grpID_to_eleID, %grpLen); 
	while (my ($t_EID, $t_GID) = each %eleID_to_grpID) {
		push(@{$grpID_to_eleID{$t_GID}}, $t_EID); 
		$grpLen{$t_GID} ++; 
	}
	for my $t_GID (sort { $grpLen{$b}<=>$grpLen{$a} || $a <=> $b } keys %grpID_to_eleID) {
		push(@back_array, [ @{$grpID_to_eleID{$t_GID}} ]); 
	}
	
	return (\@back_array); 
}# divide_group_fromArray() 


=head1 divide_group( \%relation_pairs )

Required   : Do be aware that in %relation_pairs, each pair should have two formed key pairs, {k1}{k2} and {k2}{k1}

Input      : Format of %relation_pairs is 
             {$element_1}{$element_2} = 1; 
             {$element_2}{$element_1} = 1; 
             {$element_3}{$element_4} = 1; 
             {$element_4}{$element_3} = 1; 
             {$element_5}{$element_6} = 1; 
             {$element_5}{$element_5} = 1; 
             {$element_1}{$element_6} = 1; 
             {$element_6}{$element_1} = 1; 
             ...... 

Function   : Divide elements in %relation_pairs into groups. 
             Each group contains related elements, and there is no related pairs between groups. 

Output     : (\@back_groups_array)
 @back_groups_array format is ( [group_1_ele_1, group1_ele_2, ...], [group_2_ele_1, group_2_ele_2, ...], ... ); 

=cut
sub divide_group {
	my ( $pair_href ) = @_; 
	my %new_excl ; 
	my @back_array; 
	for my $tk1 (keys %$pair_href) {
		defined $new_excl{$tk1} and next; 
		my ($b_href, $e_href) = &extract_group( $pair_href, \%new_excl, $tk1 ); 
		push(@back_array, [ sort keys %$b_href ]); 
		map { $new_excl{$_} = 1; } @{$back_array[-1]}; 
	}
	return( \@back_array ); 
}# divide_group() 

=head1 extract_group( \%relation_pairs, \%excluded_IDs, $in_seed_element ) 

Input      : Format of %relation_pairs is 
             {$element_1}{$element_2} = 1; 
             {$element_3}{$element_4} = 1; 
             {$element_5}{$element_6} = 1; 
             {$element_1}{$element_6} = 1; 
             ...... 

             Format of %excluded_IDs is [IDs that should not be considered.]
             {$element_x1} = 1; 
             {$element_x2} = 1; 
             ...... 

Function   : Extract a group which contains all and only contains the elements related to $in_seed_element by some chaining. 

Output     : (\%group_IDs, \%new_excluded_IDs)
 %group_IDs is in format {$in_seed_element} = 1, {$relate_ID_1} = 1, {$relate_ID_2} = 1, ...
 %new_excluded_IDs is in the same format of %exluded_IDs, and include all IDs in %excluded_IDs and IDs in %group_IDs; 

=cut
sub extract_group {
	my ( $pair_href, $excl_href, $in_key ) = @_; 
	my %back_hash; 
	my %new_excl = %$excl_href; 
	$back_hash{$in_key} = 1; 
	$new_excl{$in_key} = 1; 
	for my $tk ( keys %{$pair_href->{$in_key}} ) {
		defined $new_excl{$tk} and next; 
		$back_hash{$tk} = 1; 
		$new_excl{$tk} = 1; 
		my ($b_href, $e_href) = &extract_group( $pair_href, \%new_excl, $tk); 
		map { $new_excl{$_} = 1; } keys %$e_href; 
		map { $back_hash{$_} = 1; } keys %$b_href; 
	}
	return (\%back_hash, \%new_excl); 
}# extract_group() 

=head1 extract_group_fromHash( \%relation_pairs, \%excluded_IDs, $in_seed_element ) 

Input      : Format of %relation_pairs is 
             {$element_1}{$element_2} = 1; 
             {$element_3}{$element_4} = 1; 
             {$element_5}{$element_6} = 1; 
             {$element_1}{$element_6} = 1; 
             ...... 

             Format of %excluded_IDs is [IDs that should not be considered.]
             {$element_x1} = 1; 
             {$element_x2} = 1; 
             ...... 

Function   : Extract a group which contains all and only contains the elements related to $in_seed_element by some chaining. 

Output     : (\%group_IDs, \%new_excluded_IDs)
 %group_IDs is in format {$in_seed_element} = 1, {$relate_ID_1} = 1, {$relate_ID_2} = 1, ...
 %new_excluded_IDs is in the same format of %exluded_IDs, and include all IDs in %excluded_IDs and IDs in %group_IDs; 

=cut
sub extract_group_fromHash {
	my ( $pair_href, $excl_href, $in_key ) = @_; 
	$excl_href //= {}; 
	defined $in_key or die "Must give \$in_seed_element for mathSunhh::extract_group_fromHash()\n"; 
	my (%back_hash, %new_excl); 
	%new_excl = %$excl_href; 

	my @array_input; 
	for my $tk1 (keys %$pair_href) {
		defined $new_excl{$tk1} and next; 
		for my $tk2 ( keys %{$pair_href->{$tk1}} ) {
			defined $new_excl{$tk2} and next; 
			push(@array_input, [ $tk1, $tk2 ]); 
		}
	}
	my ($back_href, $new_excl_href) = &extract_group_fromArray( \@array_input, \%new_excl, $in_key ); 
	
	return($back_href, $new_excl_href); 
}# extract_group_fromHash() 

=head1 extract_group_fromArray( \@elements_in_groups, \%excluded_IDs, $in_seed_element ) 

Input      : Format of @elements_in_groups is 
             ([ ele1, ele2, ... ], 
             [ele3, ele4, ele5, ...], 
             [ele1, ele3], 
             [ele6, ele7, ele8], 
             [ele11, ele10, ele9], 
             ......
             )

             Format of %excluded_IDs is [IDs that should not be considered.]
             {$element_x1} = 1; 
             {$element_x2} = 1; 
             ...... 

Function   : Extract a group which contains all and only contains the elements related to $in_seed_element by some chaining. 

Output     : (\%group_IDs, \%new_excluded_IDs)
 %group_IDs is in format {$in_seed_element} = 1, {$relate_ID_1} = 1, {$relate_ID_2} = 1, ...
 %new_excluded_IDs is in the same format of %exluded_IDs, and include all IDs in %excluded_IDs and IDs in %group_IDs; 

=cut
sub extract_group_fromArray {
	my ( $pair_aref, $excl_href, $in_key ) = @_; 
	$excl_href //= {}; 
	defined $in_key or die "Must give \$in_seed_element for mathSunhh::extract_group_fromArray()\n"; 
	my (%back_hash, %new_excl); 
	%new_excl = %$excl_href; 

	my @array_input; 
	for my $ar1 (@$pair_aref) {
		my @ta; 
		for my $eleID (@$ar1) {
			defined $new_excl{$eleID} and next; 
			push(@ta, $eleID); 
		}
		scalar(@ta) > 0 and push(@array_input, [@ta]); 
	}
	my ($back_aref) = &divide_group_fromArray(\@array_input); 
	for my $ar1 (@$back_aref) {
		my $is_has = 0; 
		for my $eleID (@$ar1) {
			$eleID eq $in_key or next; 
			$is_has = 1; 
			last; 
		}
		if ( $is_has == 1 ) {
			for my $eleID (@$ar1) {
				$back_hash{$eleID} = 1; 
				$new_excl{$eleID} = 1; 
			}
			last; 
		}
	}
	
	return(\%back_hash, \%new_excl); 
}# extract_group_fromArray() 

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
	my $loci_2 = &mergeLocBlk('', $aref_2, 'dist2join' => 1); 
	for my $ar2 (@$loci_2) {
		my ( $t_specSE, $t_ovlpSE, $t_idx_start ) = &sep_loc2_by_loc1_singleLoc2( $aref_1, $ar2, $aref_colN ); 
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

Function   : According to \@reference_1to2, map loc_1 positions [ $targetS_1, $targetE_1 ] to loc_2 positions, and return them. 

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

=head1 minmax(\@numbers)

Return      : ($min, $max)

=cut
sub minmax {
	my ($back_min, $back_max); 
	for my $t (@{$_[0]}) {
		$back_min //= $t; 
		$back_max //= $t; 
		$back_min > $t and $back_min = $t; 
		$back_max < $t and $back_max = $t; 
	}
	return($back_min, $back_max); 
}# minmax 

=head1 log10( $number )

Return       : ( $number_value )

=cut
sub log10 {
	my $n = shift; 
	return(log($n)/log(10)); 
}# log10() 

############################################################
#  Sub-routines for number and loci indexing. 
############################################################
=head1 _decimal_to_hexa ( $number_in_decimal )

Please don't provide decimal fraction!!! 

Return        : ( $number_in_hexadecimal )
=cut
sub _decimal_to_hexa {
	my $add = ''; 
	$_[0] < 0 and do { $add = '-'; $_[0] = abs($_[0]); }; 
	return( $add . sprintf("%X", $_[0]) ); 
} # _decimal_to_hexa () 

=head1 _hexa_to_decimal ( $number_in_hexa )

Return        : ( $number_in_decimal )
=cut
sub _hexa_to_decimal {
	my $add = ''; 
	if ( $_[0] =~ s!^\-!! ) {
		$add = '-'; 
	}
	$_[0] =~ /\A(?:0?[xX])?(?:_?[0-9a-fA-F])*\z/ or &stopErr("[Err] Input [$_[0]] is not a valid hex digit string.\n"); 
	return( $add . hex($_[0]) ); 
} # _hexa_to_decimal () 

=head1 _Encode_deciN ( $Number, $max_position_len'Default:11' )

Return        : ( \@num_in_each_position )

Description   :

                I am still used to this forward direction instead of reverse direction. 
                Divide input $Number 13456 with $max_position_len==11 into [0,0,0,0,0,0,1,3,4,5,6].
                $Number must be a positive integer.
                Output of _Encode_deciN( 8123456789012 , 11 ) will be [qw/812 3 4 5 6 7 8 9 0 1 2/]

=cut
sub _Encode_deciN {
        my ($n, $maxP) = @_;
        $maxP //= 11;
        my @back = (0) x $maxP;
        my $i = -1;
        my $t_n;
        for (my $i=$maxP-1; $i >= 0 and $n > 0; $i--) {
                $t_n = int($n/10);
                $back[$i] = $n - $t_n*10;
                $n = $t_n;
        }
	if ( $n > 0 ) {
		$back[0] += $n*10; 
	}

        return(\@back);
}# _Encode_deciN ()

=head1 _Decode_deciN ( \@deciN_array )

Return        : ( $decimal_integer )

Description   : Convert [qw/2 1 0 9 8 7 6 5 4 3 2 810] to number 8123456789012 . 

=cut
sub _Decode_deciN {
	my ($ar) = @_; 
	my $back = 0; 
	my $i=-1; 
	for my $t (reverse @$ar) {
		$i++; 
		$back += 10 ** $i * $t; 
	}

	return($back); 
}# _Decode_deciN ()

=head1 index_SEloc( [ [$loc_1_S, $loc_1_E], [$loc_2_S, $loc_2_E], ... ], {'binSize' => 1000 , 'maxP' => 11 } )

Return        : (\%db_loc)

Usage case 1  : 
    my @loc = ( [1,100], [2 , 384745], [10,21], [50,60], [384745, 984672], [20,55], [1, 9999], [4,11], [10,20], [10,21], [10,19]); 
    my %locDb_h = %{ &index_SEloc( \@loc, { 'binSize' => 100 , 'maxP' => 11 } ) }; 
    print "intact raw data Start:\n"; 
    for (my $i=0; $i<@loc; $i++) {
      print join("\t", "$loc[$i][0]-$loc[$i][1]", "$locDb_h{'realLoc'}[$i][0]-$locDb_h{'realLoc'}[$i][1]")."\n"; 
    }
    print "intact raw data End:\n"; 


    ## Check if it is well sorted : 
    print "newly sorted Start :\n"; 
    my @sorted_realLoc = @{ &ret_sortedRealLoc( \%locDb_h ) }; 
    for my $ar1 (@sorted_realLoc) {
      print join("\t", @$ar1)."\n"; 
    }
    print "newly sorted End :\n"; 

    ## Find overlapping realLocs with given SEloc
    my @seLoc = (21, 30); 
    @ARGV >= 2 and @seLoc = @ARGV[0,1]; 
    print "\nRegion to overlap to: $seLoc[0] - $seLoc[1]\n\n"; 
    my @loc_realIdx = @{ &_map_loc_to_realIdx( \%locDb_h, \@seLoc ) }; 
    my @loc_realLoc = @{$locDb_h{'realLoc'}}[@loc_realIdx]; 
    for (my $i=0; $i<@loc_realIdx; $i++) {
      print join("\t", "realIdx:$loc_realIdx[$i]", "realLoc:_$loc_realLoc[$i][0]_-_$loc_realLoc[$i][1]")."\n";
    }



Description   : This is the main function I want to use in RNAseq counting. 
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
	($back{'minV'}, $back{'maxV'}) = &minmax( [ map { ($_->[0], $_->[1]) } @{$back{'realLoc'}} ] ); 
	&_binRealLoc2BinLoc( \%back ); 
	&_fill_rawN2rawIdx( \%back ); 
	&_fill_struct_SEloc( \%back ); 
	&_fill_sorted_rawIdx( \%back ); 
	
	return(\%back); 
}# index_SEloc ()

=head1 _map_loc_to_realIdx ( \%db_loc, [$chkLoc_S, $chkLoc_E] ) 

Return       : ( [ $overlap_realLocIdx_1, $overlap_realLocIdx_2, ... ] )

Description  : $overlap_realLocIdx can be used in $db_loc{'realLoc'}[$overlap_realLocIdx] (Value=[$realS, $realE]); 

=cut
sub _map_loc_to_realIdx {
	my ($locDB_hr, $locSE_ar) = @_; 
	my @back; 
	my @loc_rawIdx = @{ &_map_loc_to_rawIdx( $locDB_hr , $locSE_ar ) }; 
	my ($tS,$tE) = @$locSE_ar; 
	my @realIdx_arr = @{ &_locDb_rawIdx2realIdx( $locDB_hr, \@loc_rawIdx, 1 ) }; 
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
	my @sorted_realIdx = @{ &_locDb_rawIdx2realIdx( $ah, \@sorted_rawIdx, 1 ) }; 
	@back = map { [@$_] } @{$ah->{'realLoc'}}[ @sorted_realIdx ]; 

	return(\@back); 
}# ret_sortedRealLoc () 

=head1 index_numbers( \@raw_numbers, { 'maxP' => 11 } )

Return       : (\%db_num)

Function     : Index the numbers, producing keys in \%db_num, including : {'maxP|rawNum|rawN2rawIdx|struct|sorted_rawIdx|rawN2srtIdx'}
    ->{'maxP'}          =>  11 , or someother number (>0) showing the depth of structure and the depth (length)  of encoded number path (for &_Encode_deciN())
    ->{'rawNum'}        =>  [@raw_numbers]
    ->{'rawN2rawIdx'}   =>  { $rawNum => [ $rawIdx1, $rawIdx2, ... ] } according to ->{'rawNum'}
    ->{'struct'}        =>  An array of arrays with depth (rows) of ->{'mapP'} and width (cols) of [0-9] or the maximum of the first(0) row. 
    ->{'sorted_rawIdx'} =>  [ ->{'rawN2rawIdx'}{$rawNum_1}, {'rawN2rawIdx'}{$rawNum_2}, ... ] 
                            === [ [$rawIdx3, $rawIdx8, ...], [$rawIdx2], ... ], 
                               # Here, $rawIdx3_rawNumber ( === ->{'rawNum'}[$rawIdx3] ) is equal to $rawIdx8_rawNumber, and smaller than $rawIdx2_rawNumber;
    ->{'rawN2srtIdx'}   =>  { $rawNum => $index_of_'sorted_rawIdx'_array }

Usage case 1 : 
    my @numbers = ( 1.. 100, 40 .. 50 ,105, 103, 105, 150, 1, 999, 3000, 888, 444, 999, 3000, 999, 99999999999999, 111, 99999999999399 ); 
    my %numDb_h = %{ &index_numbers( \@numbers, { 'maxP' => 11 } ) }; 

    print "input\@=@numbers\n"; 
    print "stored=@{$numDb_h{'rawNum'}}\n"; 
    ## Get sorted numbers : 
    print "sorted=" . join(" ", @{ &ret_sortedNum( \%numDb_h )}) . "\n"; 
    ## Get sorted raw indices : 
    my @sorted_rawIdx = map { @$_ } @{$numDb_h{'sorted_rawIdx'}}; 
    for my $rawIdx ( @sorted_rawIdx ) {
      print join("\t", "rawIdx:$rawIdx", "rawNum:$numDb_h{'rawNum'}[$rawIdx]")."\n"; 
    }

    ## Locate a number in database : 
    my $number_to_check = 99999999999300; 
    print "searching=$number_to_check\n"; 
    my $srtIdx = &_map_number_to_srtIdx( \%numDb_h, $number_to_check ); 
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
	&_fill_rawN2rawIdx( \%back ); 
	&_fill_struct( \%back ); 
	&_fill_sorted_rawIdx( \%back ); 

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
				return( &_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $i-1], $j ], 'max', $path_depth ) ); 
			}
			# S2: There is no idxV smaller than input in current level ($i), so find the biggest in the most recent uppper level. 
			#     I don't want to combine this and the S1 step, because I want to use &_defined_pathInStruct() to confirm there is no problem. 
			for (my $lvl=$i-1; $lvl >= 0; $lvl--) {
				my $lvl_struct = &_defined_pathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl] ] ) or &stopErr("[Err] Weired get here! 'left'\n"); 
				my $lvl_idx    = $ar_path_in->[$lvl]; 
				# Check if this level ($lvl) has an idxV < $ar_path_in->[$lvl]; 
				for ( my $j=$lvl_idx-1; $j>=0; $j-- ) {
					defined $lvl_struct->[$j] or next; 
					return( &_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl], $j ], 'max', $path_depth ) ); 
				}
			}
			# S3: There is no smaller path found in total dataset, so I should use the smallest in database. 
			return( &_get_pathInStruct( $ar_struct, 'min', $path_depth ) ); 
		} elsif ( $side eq 'right' ) {
			# S1: Check if the current level ($i) has a idxV > $ar_path_in->[$i]; 
			for (my $j=$ar_path_in->[$i]+1; $j<=@$tmp_ar; $j++) {
				defined $tmp_ar->[$j] or next; 
				# If has such an element
				return( &_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $i-1], $j ], 'min', $path_depth ) ); 
			}
			# S2: Continue to check from the most recent upper level ($i-1). 
			for (my $lvl=$i-1; $lvl >= 0; $lvl--) {
				my $lvl_struct = &_defined_pathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl] ] ) or &stopErr("[Err] Wiered get here! 'right'\n"); 
				my $lvl_idx    = $ar_path_in->[$lvl]; 
				# Check if this level contains an idxV > $ar_path_in->[$lvl]; 
				for ( my $j=$lvl_idx+1; $j<=@$lvl_struct; $j++ ) {
					defined $lvl_struct->[$j] or next; 
					return( &_get_subPathInStruct( $ar_struct, [ @{$ar_path_in}[0 .. $lvl], $j ], 'min', $path_depth ) ); 
				}
			}
			# S3: There is no smaller path found in total dataset, so I should use the biggest in database. 
			return( &_get_pathInStruct( $ar_struct, 'max', $path_depth ) ); 
		} elsif ( $side eq 'both' ) {
			# S1: Find the left most path 
			my $path_left = &_sideClose_pathInStruct( $ar_struct, $ar_path_in, 'left' ); 
			my $dist_left = abs( &_distance_paths( $ar_path_in, $path_left ) ); 
			# S2: Find the right most path 
			my $path_right = &_sideClose_pathInStruct( $ar_struct, $ar_path_in, 'right' ); 
			my $dist_right = abs( &_distance_paths( $path_right, $ar_path_in ) ); 
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
	my $back = &_Decode_deciN($p1) - &_Decode_deciN($p2); 
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
			my $sub_back = &_lsAll_inStruct( $ar_struct->[$i], $maxP-1 ); 
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
	my $num_path = &_Encode_deciN($number); 
	my $closest_path = &_sideClose_pathInStruct( $ah->{'struct'}, $num_path, $side ); 
	my $closest_rawN = $ah->{'rawNum'}[ &_get_valueFromStruct( $ah->{'struct'}, $closest_path )->[0] ]; 

	return( $closest_rawN ); 
}# _map_number_to_rawNum ()

=head1 _map_number_to_srtIdx( \%db_number, $number_to_locate, 'both|left|right' )

Return       : ($idx_in_'sorted_rawIdx_array')

Description  : The third parameter will be 'both' in default. 

=cut
sub _map_number_to_srtIdx {
	my ($ah, $number, $side) = @_; 
	$side //= 'both'; 
	my $closest_rawN = &_map_number_to_rawNum( $ah, $number, $side ); 
	# my $closest_rawN = &_Decode_deciN($closest_path) ; 
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

	my $sub_struct = &_defined_pathInStruct( $ar_struct, $root_path ) or &stopErr("[Err] No structure exists\n"); 
	my $sub_path   = &_get_pathInStruct( $sub_struct, $mm, $path_depth-scalar(@$root_path) ); 
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
                @encoded_num comes from &_Encode_deciN() function, and it is treated as a full path in @struct. 
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

Description   : Change 'rawIdx' in database to 'realIdx' numbers; 
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

FUnction      : 
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

Function    : Add ->{'struct'}     => [ structure_database ]
              final value of 'struct' is ->{'rawN2rawIdx'}{$rawN} 

=cut
sub _fill_struct {
	my ( $ah ) = @_; 
	$ah->{'struct'} = []; 
	for my $rawN (keys %{$ah->{'rawN2rawIdx'}}) {
		&_put_num2Struct( $ah->{'struct'}, &_Encode_deciN($rawN, $ah->{'maxP'}), $ah->{'rawN2rawIdx'}{$rawN} ); 
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
	my %eLoc_db = %{ &index_numbers(\@eLoc, { 'maxP' => $ah->{'maxP'} }) }; 
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
	my %sLoc_db = %{ &index_numbers(\@sLoc, { 'maxP' => $ah->{'maxP'} } ) }; 
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
	for my $ar_rawIdx_INsLocDb ( map { $_->[1] } @{ &_lsAll_inStruct( $ah->{'struct'} , $ah->{'maxP'} ) } ) {
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
	for my $ar ( @{ &_lsAll_inStruct( $ah->{'struct'}, $ah->{'maxP'} ) } ) {
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
	my $rawNum_toChk_ar = $ms_obj->map_windows( 'position' => $posi, 'wind_hash' => $ah->{'wind'} ); 
	my @back; 
	my %has; 
	for my $rawNum_toChk ( @$rawNum_toChk_ar ) {
		defined $ah->{'rawN2rawIdx'}{$rawNum_toChk} or next; 
		# my $realIdx_ar = &_locDb_rawIdx2realIdx( $ah, $ah->{'rawN2rawIdx'}{$rawNum_toChk}, 1 ); 
		for my $t ( @{ &_locDb_rawIdx2realIdx( $ah, $ah->{'rawN2rawIdx'}{$rawNum_toChk}, 1 ) } ) {
			defined $has{ $t } and next; 
			$has{ $t } = 1; 
			$posi >= $ah->{'realLoc'}[$t][0] and $posi <= $ah->{'realLoc'}[$t][1] and push(@back, $t); 
		}
	}

	return(\@back); 
}# _hasPos_inLocDb ()


=head1 _map_loc_to_rawIdx ( \%db_loc, [$chkLoc_S, $chkLoc_E] )

Return        : ([$idx1_in_'rawLoc_array', $idx2_in_'rawLoc_array', ...]) 

Description   : Better ignore this inner function. 
                $idx1_in_'rawLoc_array' can be used in $db_loc{'rawLoc'}[$idx1_in_rawLoc_array] (Value=[$rawLocS, $rawLocE]). 
                Be aware that when there are overlaps between rawLocs, the return might be not complete!!! 

=cut
sub _map_loc_to_rawIdx {
	my ($locDB_hr, $locSE_ar) = @_; 
	my @locSE = (@$locSE_ar); 
	$locSE[0] > $locSE[1] and @locSE = @locSE[1,0]; 
	my @back; 

	my $find_idxS = &_map_number_to_srtIdx( $locDB_hr, $locSE[0], 'both' ); 
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
    Add ->{'wind'}                => $ms_obj->setup_windows(), in which ->{'wind'}{'loci'}
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
	my %locDb = %{ &index_SEloc( \@realNum, { 'binSize' => 0 } ) }; 

	# Construct windows (bins)
	my %wind = %{ $ms_obj->setup_windows( 'ttl_start'=>$ah->{'minV'}, 'ttl_end'=>$ah->{'maxV'}, 'wind_size'=>$ah->{'binSize'}, 'wind_step'=>$ah->{'binSize'}, 'minRatio'=>0 ) }; 
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
			my @idxS_arr = @{ $ms_obj->map_windows( 'position' => $real_SE[0], 'wind_hash'=>\%wind ) }; 
			my @idxE_arr = @{ $ms_obj->map_windows( 'position' => $real_SE[1], 'wind_hash'=>\%wind ) }; 
			my ($si_S, $si_E) = &minmax( [@idxS_arr, @idxE_arr] ); 
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

