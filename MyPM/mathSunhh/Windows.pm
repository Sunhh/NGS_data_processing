package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 
use Scalar::Util; 


=head2 setup_windows( 'ttl_start'=>1, 'ttl_end'=>99999999, 'wind_size'=>1000, 'wind_step'=> // 'wind_size', 
  'minRatio'=> 0, 
  'max_end_wind_num' => 0
)

Required: 
 

Function: return a list of windows in hash reference according to given [window_size, window_step, total_start, total_end]

  If 'max_end_wind_num' is 0, any window accepted by minRatio are kept. 
  If 'max_end_wind_num' > 0, maximum 'max_end_wind_num' number of window reaching the end position are kept. 

Return  : \%back_wind; 

 'keys'  => values
 'info' => {'ttl_start/ttl_end/wind_size/wind_step/minRatio/windSloci' => values} 
   {'info'}{'windSloci'} = [ start_pos_0, start_pos_1, start_pos_2, ... ]
   {'info'}{'Sloc2wi'}{$start_pos} => window_idx_in_windSloci 
 'loci'  => {
              'wind_start_position' => [start_pos, end_pos, interval_len]
            }
 
=cut
sub setup_windows {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	$parm{'ttl_start'} = $parm{'ttl_start'} // 1; 
	$parm{'ttl_end'}   = $parm{'ttl_end'}   // 99999999; 
	$parm{'wind_size'} = $parm{'wind_size'} // 1000; 
	$parm{'wind_step'} = $parm{'wind_step'} // $parm{'wind_size'}; 
	$parm{'minRatio'}  = $parm{'minRatio'}  // 0; 

	$parm{'max_end_wind_num'}= $parm{'max_end_wind_num'} // 0; 
	my $end_num = 0; 

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
		$ei == $parm{'ttl_end'} and $end_num ++; 
		my $cur_len = $ei-$si+1; 
		$back_wind{'loci'}{$si} = [$si, $ei, $cur_len]; 
		push(@{$back_wind{'info'}{'windSloci'}}, $si); 
		$back_wind{'info'}{'Sloc2wi'}{$si} = $#{$back_wind{'info'}{'windSloci'}}; 

		$parm{'max_end_wind_num'} > 0 and $end_num >= $parm{'max_end_wind_num'} and last; 
	}
	return \%back_wind; 
}# sub setup_windows

=head2 map_windows ( 'posi/position'=>Integer, 'wind_hash'=>setup_windows->(), 
  'ttl_start'=>1, 'ttl_end'=>99999999, 'wind_size'=>1000, 'wind_step'=>'wind_size', 'minRatio'=>0, 
  'max_end_wind_num' => 0 
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
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	my $posi = $parm{'posi'} // $parm{'position'} // &stopErr("[Err] No position assigned.\n"); 
	$parm{'max_end_wind_num'}= $parm{'max_end_wind_num'} // 0; 
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
	my $end_num = 0; 
	@back_si = reverse(@back_si); 
	my @t1; 
	for my $ti (@back_si) {
		$ti+$parm{'wind_size'}-1 >= $parm{'ttl_end'} and $end_num ++; 
		push(@t1, $ti); 
		$parm{'max_end_wind_num'} > 0 and $end_num >= $parm{'max_end_wind_num'} and last; 
	}
	@back_si = @t1; 

	
	return \@back_si; 
}# sub map_windows 

1; 
