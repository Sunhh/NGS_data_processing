package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 



=head1 complete_pair_for_grouping ( \%relation_pairs )

Input      : 

 Format of %relation_pairs is 
   {$element_1}{$element_2} = 1; 
   {$element_3}{$element_4} = 1; 
   {$element_4}{$element_3} = 1; 
   {$element_5}{$element_6} = 1; 
   {$element_1}{$element_6} = 1; 
   {$element_6}{$element_1} = 1; 
   ...... 

Output     : 

 ( \%complete_relation_pairs )
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

Input      : 

 Format of %relation_pairs is 
   {$element_1}{$element_2} = 1; 
   {$element_2}{$element_1} = 1; 
   {$element_3}{$element_4} = 1; 
   {$element_4}{$element_3} = 1; 
   {$element_5}{$element_6} = 1; 
   {$element_5}{$element_5} = 1; 
   {$element_1}{$element_6} = 1; 
   {$element_6}{$element_1} = 1; 
   ...... 

Function   : 

 Divide elements in %relation_pairs into groups. 
   Each group contains related elements, and there is no related pairs between groups. 
   This should be much faster than &mathSunhh::divide_group(); 

Output     : (\@back_groups_array)
 
 @back_groups_array format is ( [group_1_ele_1, group1_ele_2, ...], [group_2_ele_1, group_2_ele_2, ...], ... ); 

=cut
sub divide_group_fromHash {
	my ( $pair_href ) = @_; 
	my @array_input; 
	map { push( @array_input, [ $_, keys %{$pair_href->{$_}} ] ); } keys %$pair_href; 
	my ($back_aref) = &mathSunhh::divide_group_fromArray(\@array_input); 
	
	return ($back_aref); 
}# divide_group_fromHash() 

=head1 divide_group_fromArray ( \@elements_in_groups ) 

Required   : 

Input      : 

 Format of @elements_in_groups is 
   ([ ele1, ele2, ... ], 
    [ele3, ele4, ele5, ...], 
    [ele1, ele3], 
    [ele6, ele7, ele8], 
    [ele11, ele10, ele9], 
    ......)

Function   :

 Divide elements in @elements_in_groups into groups. 
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

Input      : 

 Format of %relation_pairs is 
   {$element_1}{$element_2} = 1; 
   {$element_2}{$element_1} = 1; 
   {$element_3}{$element_4} = 1; 
   {$element_4}{$element_3} = 1; 
   {$element_5}{$element_6} = 1; 
   {$element_5}{$element_5} = 1; 
   {$element_1}{$element_6} = 1; 
   {$element_6}{$element_1} = 1; 
   ...... 

Function   : 

 Divide elements in %relation_pairs into groups. 
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
		my ($b_href, $e_href) = &mathSunhh::extract_group( $pair_href, \%new_excl, $tk1 ); 
		push(@back_array, [ sort keys %$b_href ]); 
		map { $new_excl{$_} = 1; } @{$back_array[-1]}; 
	}
	return( \@back_array ); 
}# divide_group() 

=head1 extract_group( \%relation_pairs, \%excluded_IDs, $in_seed_element ) 

Input      :

 Format of %relation_pairs is 
   {$element_1}{$element_2} = 1; 
   {$element_3}{$element_4} = 1; 
   {$element_5}{$element_6} = 1; 
   {$element_1}{$element_6} = 1; 
   ...... 

 Format of %excluded_IDs is [IDs that should not be considered.]
   {$element_x1} = 1; 
   {$element_x2} = 1; 
   ...... 

Function   : 

 Extract a group which contains all and only contains the elements related to $in_seed_element by some chaining. 

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
		my ($b_href, $e_href) = &mathSunhh::extract_group( $pair_href, \%new_excl, $tk); 
		map { $new_excl{$_} = 1; } keys %$e_href; 
		map { $back_hash{$_} = 1; } keys %$b_href; 
	}
	return (\%back_hash, \%new_excl); 
}# extract_group() 

=head1 extract_group_fromHash( \%relation_pairs, \%excluded_IDs, $in_seed_element ) 

Input      : 

 Format of %relation_pairs is 
   {$element_1}{$element_2} = 1; 
   {$element_3}{$element_4} = 1; 
   {$element_5}{$element_6} = 1; 
   {$element_1}{$element_6} = 1; 
   ...... 

 Format of %excluded_IDs is [IDs that should not be considered.]
   {$element_x1} = 1; 
   {$element_x2} = 1; 
   ...... 

Function   : 

 Extract a group which contains all and only contains the elements related to $in_seed_element by some chaining. 

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
	my ($back_href, $new_excl_href) = &mathSunhh::extract_group_fromArray( \@array_input, \%new_excl, $in_key ); 
	
	return($back_href, $new_excl_href); 
}# extract_group_fromHash() 

=head1 extract_group_fromArray( \@elements_in_groups, \%excluded_IDs, $in_seed_element ) 

Input      : 

  Format of @elements_in_groups is 
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

Function   :

 Extract a group which contains all and only contains the elements related to $in_seed_element by some chaining. 

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
	my ($back_aref) = &mathSunhh::divide_group_fromArray(\@array_input); 
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

1; 
