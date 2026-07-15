package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 



=head2 repArr (\@units, 'each'=>1, 'times'=>1, 'length'=>\@units) 

Required : \@units

Function : 

	A repeat function similar to rep() in R. 
	First apply 'each', then apply 'times', and at least apply 'length'; 

Return   : \@repeated_eles 

=cut
sub repArr {
	my $ref_arr = shift; 
	my %parm = &mathSunhh::_setHashFromArr(@_); 
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
		for ( &mathSunhh::permutations(\@rest, $n-1) ) {
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
		for (&mathSunhh::combinations(\@rest, $n-1)) {
			push(@comb, [$val, @$_]); 
		}
	}
	return @comb; 
}#sub combinations

=head1 dvd_array( \@array_to_be_divided, $number_of_subgroups, $If_keep_order[default=0], $prev_grp_colN )

Function : Divide @array_to_be_divided into $number_of_subgroups subgroups. 

 If $If_keep_order is TRUE, the elements in subgroups will be sequenctial same to the raw order. 
 If $If_keep_order is TRUE and $prev_grp_colN is not N, $prev_grp_colN should be a colNumber and each sub-group has only one type of char in $array_to_be_divided[$prev_grp_colN];


Return   : ( [ \@subgroup_1, \@subgroup_2, ... ] )

Example  : 
  &mathSunhh::dvd_array( [0,1,2,3,4,5,6,7] , 3 )    returns [ [0,3,6], [1,4,7], [2,5] ] ; 
  &mathSunhh::dvd_array( [0,1,2,3,4,5,6,7] , 3, 1 ) returns [ [0,1,2], [3,4,5], [6,7] ] ; 

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

=head1 randSlct_num( $total_num, $slct_Num )

Return      : ( \@slct_indice_sorted )

Annotation  : There is no replacement in the randomly selecting; 

=cut
sub randSlct_num {
	my ($total_num, $n) = @_; 
	my @back; 
	$n = int($n); 
	$n > 0 or return(\@back); 
	my @idx = ( 0 .. ($total_num-1) ); 
	$total_num <= $n and return(\@idx); 
	while ( @idx > 0 ) {
		my $j = rand( $#idx + 1 ); 
		push(@back, splice( @idx, $j, 1 )); 
		$#back+1 >= $n and last; 
	}
	@back = sort { $a <=> $b } @back; 
	return(\@back); 
}# sub randSlct_num () 

=head1 create_randNum ( $num_digits )

Return      : ( $a_number )

=cut
sub create_randNum {
	my $num_digits = shift; 
	my $nn = ""; 
	for (1 .. $num_digits) {
		my $r=int(rand(10));
		$r == 10 and $r = 9;
		$nn .= $r;
	}
	$nn =~ s!^0!1!;
	return($nn);
}# sub create_randNum() 

1; 
