package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 


############################################################
#  Sub-routines for number and loci indexing. 
############################################################
=head1 _decimal_to_hexa ( $number_in_decimal )

 Please don't provide decimal fraction!!! 

Return        : ( $number_in_hexadecimal )

  _decimal_to_hexa( -100 ) returns ( '-64' )
  _decimal_to_hexa( 200, -100 ) returns ( 'C8', '-64' )

=cut
sub _decimal_to_hexa {
	my @back; 
	scalar(@_) == 0 and return; 
	for my $t1 (@_) {
		my $t = $t1; 
		my $add = ''; 
		$t < 0 and do { $add = '-'; $t = abs($t); }; 
		push( @back, $add . sprintf("%X", $t) ); 
	}
	if (scalar(@_) > 1) {
		return (@back); 
	} else {
		return ($back[0]); 
	}
} # _decimal_to_hexa () 

=head1 _hexa_to_decimal ( $number_in_hexa )

Return        : ( $number_in_decimal )

  _hexa_to_decimal( -64 ) returns ( -100 )
  _hexa_to_decimal( -64, 'C8') returns ( -100, 200 )

=cut
sub _hexa_to_decimal {
	my @back; 
	scalar(@_) == 0 and return; 
	for my $t1 (@_) {
		my $t = $t1; 
		my $add = ''; 
		$t =~ s!^\-!! and $add = '-'; 
		$t =~ /\A(?:0?[xX])?(?:_?[0-9a-fA-F])*\z/ or &stopErr("[Err] Input [$t] is not a valid hex digit string.\n"); 
		push(@back, $add . hex($t)); 
	}
	if (scalar(@_) > 1) {
		return (@back); 
	} else {
		return ($back[0]); 
	}
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

1; 
