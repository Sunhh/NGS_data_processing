package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 


############################################################
#  Sub-routines. 
############################################################


=head1 get_xy_byScale( 'start_xy'=>[$sX,$sY], 'end_xy'=>[$eX,$eY], 'se_position'=>[$sPos,$ePos], 'need_position'=>[$p1,$p2,...] )

Return : ( [$needP_1_x, $needP_1_y], [$needP_2_x, $needP_2_y], ... )

=cut
sub get_xy_byScale {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
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
	$val < $min and $val = $min ; 
	$val > $max and $val = $max ; 
	my $col_i = $#$col;
	my $v = ($val-$min)/($max-$min) * $col_i;
	my $i1 = &mathSunhh::min( int($v), $col_i );
	my $i2 = &mathSunhh::min( int($v)+1, $col_i);
	my $f = $v - $i1;
	my @back;
	for (my $j=0; $j<3; $j++) {
		$back[$j] = int( $col->[$i1][$j] + $f * ($col->[$i2][$j]-$col->[$i1][$j]) );
	}
	return( sprintf("rgb(%.0f,%.0f,%.0f)", $back[0], $back[1], $back[2])  );
}# cnvt_to_rgb () 



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
			&stopErr("[Err] Unparsable column tag [$tc] in parseCol().\n"); 
		}
	}
	return (@ncols); 
}# parseCol() 
*_parseCol = \&parseCol; # backward-compatible alias for the old private name

1; 
