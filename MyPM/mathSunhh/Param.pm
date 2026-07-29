package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
# 260729: Keep empty value as 0 in parseCol; This function is moved from Color.pm;
use strict; 
use warnings; 
use LogInforSunhh; 


=head1 _setHashFromArr(@keyVal_array)
=cut

=head2 _setHashFromArr(@keyVal_array)

Required: @keyVal_array

Function: 

  @keyVal_array contain ( key1, val1, key2, val2, ... ) pairs. 
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


=head1 _addHash( 'toH'=>\%hash_ref_receptor, 'fromH'=>\%hash_ref_donor , 'replaceExist' => 0 )

Function    : 

 Add (key,value) pairs from %hash_ref_donor to %hash_ref_receptor ; 
 If 'replaceExist' is 0 , existing key in %hash_ref_receptor won't be changed, or else it will be changed. 

Return      : ( \%hash_ref_receptor )

=cut
sub _addHash {
  my %parm = &mathSunhh::_setHashFromArr(@_); 
  $parm{'toH'}   //= {}; 
  $parm{'fromH'} //= {}; 
  $parm{'replaceExist'} //= 0; 
  for my $k (keys %{$parm{'fromH'}}) {
    if ( $parm{'replaceExist'} ) {
      $parm{'toH'}{$k} = $parm{'fromH'}{$k}; 
    } else {
      $parm{'toH'}{$k} //= $parm{'fromH'}{$k}; 
    }
  }
  return($parm{'toH'}); 
}# _addHash () 

=head1 parseCol( $col_string )   (&_parseCol kept as a legacy alias)

Required    : $col_string

 $col_string : Like 0-100 | 0,3,4,5-7,1-10 | 10-8

Function    : 

 '0-5' returns (0 .. 5)
 '0,3,5-7,6-9' returns (0,3,5,6,7,6,7,8,9)
 ' -3 - -4, 10-8'  returns (-3,-4,10,9,8)

Return      : (@col_numbers)

=cut
sub parseCol {
  my @ncols; 
  for my $tc (split(/,/, shift)) {
    $tc =~ s!\s!!g; 
    # $tc =~ m!^$! and next; 
    if ( $tc =~ m!^\-?\d+$! ) {
      push(@ncols, $tc); 
    } elsif ( $tc =~ m!^(\-?\d+)\-(\-?\d+)$! ) {
      my ($s, $e) = ($1, $2); 
      if ( $s <= $e ) {
        push(@ncols, ($s .. $e)); 
      } else {
        push(@ncols, reverse($e .. $s));
      }
    } elsif ($tc eq '') {
      push(@ncols, "");
    } else {
      &stopErr("[Err] Unparsable column tag [$tc] in parseCol().\n"); 
    }
  }
  return (@ncols); 
}# parseCol() 
*_parseCol = \&parseCol; # backward-compatible alias for the old private name


1;

