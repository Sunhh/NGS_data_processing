package mathSunhh; 
# Part of mathSunhh; split out for navigability. Loaded by mathSunhh.pm (do not 'use' directly).
use strict; 
use warnings; 
use LogInforSunhh; 
use Statistics::Descriptive; 


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

=cut

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

=cut

=head2 max(@numbers)

Function: This is not a method, but a subroutine()

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
    BWA estimates the insert size distribution per 256*1024 read pairs. 
    It first collects pairs of reads with both ends mapped with a single-end quality 20 or higher and then calculates median (Q2), 
    lower and higher quartile (Q1 and Q3). 
    It estimates the mean and the variance of the insert size distribution from pairs whose insert sizes are within interval [Q1-2(Q3-Q1), Q3+2(Q3-Q1)]. 
    The maximum distance x for a pair considered to be properly paired (SAM flag 0x2) is calculated by solving equation Phi((x-mu)/sigma)=x/L*p0, 
    where mu is the mean, sigma is the standard error of the insert size distribution, L is the length of the genome, 
    p0 is prior of anomalous pair and Phi() is the standard cumulative distribution function. 
    For mapping Illumina short-insert reads to the human genome, x is about 6-7 sigma away from the mean. 
    Quartiles, mean, variance and x will be printed to the standard error output.

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
	my @clean_num;
	if (defined $r_arr and scalar(@$r_arr) > 0) {
		@clean_num = grep { !(m!^[+-]?nan$!i) } @$r_arr;
	}
	if ( (! defined $r_arr) or scalar(@clean_num) < $min_val_number) {
		for my $ta (qw/SUM COUNT MEAN MEDIAN Q1 Q3 interval_low interval_high interval_mean interval_median interval_var interval_stdev limit_low limit_high/) {
			$back{$ta} = ''; 
		}
		return \%back; 
	}
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@clean_num);
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
	for my $ta (@clean_num) {
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


=head1 p_adjust_BH( \@pvalues )  or  p_adjust_BH( \@pvalues, $n )

Function   : Benjamini-Hochberg (FDR) multiple-testing correction; a pure-Perl
             equivalent of R's p.adjust(p, method="BH", n=$n). No CPAN / R dep.

Input      : ( \@pvalues )       - raw p-values in any order.
             ( \@pvalues, $n )   - optional total count $n (mirrors R's 'n' arg;
                                    default = scalar(@pvalues); must be >= that).

Output     : ( \@adjusted )      - BH-adjusted p-values, aligned 1:1 with the
                                    input order (returned as an array reference).

Algorithm  : adj = pmin(1, cummin( n/i * p[o] )) reordered to input order, where o
             sorts p descending and i = m, m-1, ..., 1 (the ascending rank). Matches
             R's stats::p.adjust "BH" branch, including the monotone (cumulative-
             minimum) enforcement.

=cut
sub p_adjust_BH {
	my $r_p = shift; 
	ref($r_p) eq 'ARRAY' or &stopErr("[Err] p_adjust_BH() input-1st should be an array reference.\n"); 
	my $m = scalar(@$r_p); 
	$m == 0 and return []; 
	my $n = shift // $m; 
	$n >= $m or &stopErr("[Err] p_adjust_BH() n ($n) must be >= number of p-values ($m).\n"); 
	# Indices of @$r_p sorted by p-value descending (stable sort: ties keep input order).
	my @o = sort { $r_p->[$b] <=> $r_p->[$a] } 0 .. $m-1; 
	my @adj; 
	my $running_min; 
	for my $j ( 0 .. $m-1 ) {
		my $rank = $m - $j;                        # ascending rank: m, m-1, ..., 1 
		my $val  = $r_p->[$o[$j]] * $n / $rank;    # n/rank * p 
		$running_min = $val if ( !defined($running_min) or $val < $running_min ); 
		$adj[$o[$j]] = ( $running_min > 1 ) ? 1 : $running_min; 
	}
	return \@adj; 
}# p_adjust_BH() 

=head1 chisqrprob( $df, $x )

Function   : Upper-tail p-value of a chi-square distribution: P(X > $x) for $df
             degrees of freedom. Pure-Perl replacement for
             Statistics::Distributions::chisqrprob (no CPAN dep); numerically it
             equals R's pchisq($x, $df, lower.tail=FALSE).

Input      : ( $df, $x )  - degrees of freedom (>0) and the chi-square statistic.

Output     : ( $p )       - the upper-tail probability in [0,1].

Notes      : Returns 1 when $x <= 0 or $df <= 0 (conservative; matches the LRT
             guard where a non-positive statistic means "no evidence"). Computed
             as the regularized upper incomplete gamma Q($df/2, $x/2).

=cut
sub chisqrprob {
	my ($df, $x) = @_; 
	(defined $df and defined $x) or &stopErr("[Err] chisqrprob() needs ( \$df, \$x ).\n"); 
	($df <= 0 or $x <= 0) and return 1; 
	return &_gammq($df/2, $x/2); 
}# chisqrprob() 

# Regularized upper incomplete gamma  Q(a,x) = 1 - P(a,x).  (Numerical Recipes)
sub _gammq {
	my ($a, $x) = @_; 
	$x <= 0 and return 1; 
	$a >  0 or &stopErr("[Err] _gammq() needs a>0.\n"); 
	if ( $x < $a + 1 ) {
		return 1 - &_gser($a, $x);      # series gives P(a,x); Q = 1-P 
	} else {
		return &_gcf($a, $x);           # continued fraction gives Q(a,x) directly 
	}
}# _gammq() 

# Series representation of the regularized lower incomplete gamma P(a,x); x < a+1.
sub _gser {
	my ($a, $x) = @_; 
	$x <= 0 and return 0; 
	my $ITMAX = 500; 
	my $EPS   = 3e-14; 
	my $ap    = $a; 
	my $sum   = 1 / $a; 
	my $del   = $sum; 
	for ( 1 .. $ITMAX ) {
		$ap  += 1; 
		$del *= $x / $ap; 
		$sum += $del; 
		abs($del) < abs($sum) * $EPS and last; 
	}
	return $sum * exp( -$x + $a * log($x) - &_lngamma($a) ); 
}# _gser() 

# Continued-fraction representation of the regularized upper incomplete gamma Q(a,x); x >= a+1.
sub _gcf {
	my ($a, $x) = @_; 
	my $ITMAX = 500; 
	my $EPS   = 3e-14; 
	my $FPMIN = 1e-300; 
	my $b = $x + 1 - $a; 
	my $c = 1 / $FPMIN; 
	my $d = 1 / $b; 
	my $h = $d; 
	for my $i ( 1 .. $ITMAX ) {
		my $an = -$i * ($i - $a); 
		$b += 2; 
		$d = $an * $d + $b;  abs($d) < $FPMIN and $d = $FPMIN; 
		$c = $b + $an / $c;  abs($c) < $FPMIN and $c = $FPMIN; 
		$d = 1 / $d; 
		my $del = $d * $c; 
		$h *= $del; 
		abs($del - 1) < $EPS and last; 
	}
	return exp( -$x + $a * log($x) - &_lngamma($a) ) * $h; 
}# _gcf() 

# Natural log of the gamma function (Lanczos, Numerical Recipes gammln); x > 0.
sub _lngamma {
	my $xx = shift; 
	my @cof = ( 76.18009172947146, -86.50532032941677, 24.01409824083091, 
	            -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 ); 
	my $y = my $x = $xx; 
	my $tmp = $x + 5.5; 
	$tmp -= ($x + 0.5) * log($tmp); 
	my $ser = 1.000000000190015; 
	for my $j ( 0 .. 5 ) { $y += 1; $ser += $cof[$j] / $y; } 
	return -$tmp + log( 2.5066282746310005 * $ser / $x ); 
}# _lngamma() 
1; 
