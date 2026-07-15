package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 

our $_pkgObj = {}; # package-level state holder for function-form newNumber()



=head2 newNumber ( 'other_safeNumber'=>[$mathSunhhObj->{'safeNumber'}, ...], 'onlyMerge'=>0, 'debug'=>0 )

Required : Null 

Function :

	Construct 'safeNumber' hash to store not used numbers in the same object. 
	It seems that this number is not only unique in the same object, but also unique in all objects in the same module calling this mathSunhh.pm module. 

Input    : Null 

Return   : A new number not used before. 

=cut
sub newNumber {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	my $self = $_pkgObj; # function-only: the package-level number pool
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

Function : 

  Trace back all offsprings from root ID according to sub_routine_reference given. 
  Be aware that if the offspring has the same ID of rootID, this function will terminate!!! 
  This is used to avoid infinite cycles caused by relationship: rootID is a child of rootID. 

Input    : ( $rootID, sub { return @offspring; }, 'unique'=>1 )

Output   : [offspring1, offspring2, ...]

=cut
sub offspringArray {
	my $rootID = shift;  # rootID from which to find offsprings. 
	my $coderef = shift; # reference of subroutine which return a array of offspring list. 
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	$parm{'unique'} = $parm{'unique'} // 1; 

	my @cIDs;            # All Children IDs. 
	my @off1 = ( &$coderef( $rootID ) ); # Offsprings directly from rootID. 
	scalar(@off1) == 0 and return [];
	my %haveID ; 
	$haveID{$rootID} = 1; 

	for my $cID ( @off1 ) {
		defined $haveID{$cID} and next; 
		push(@cIDs, $cID);
		my @ccIDs = grep { !(defined $haveID{$_}) } @{ &mathSunhh::offspringArray( $cID, $coderef, %parm ) };
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

1; 
