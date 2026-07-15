package PopGenSunhh; 
# Self-contained population-genetics statistics (nucleotide diversity pi,
# Watterson theta, segregating sites, Tajima's D) — a drop-in replacement for
# the Bio::PopGen::Statistics functions used by reseq_tools/snpTbl_stats.pl.
#
# Data model (no BioPerl objects):
#   $inds = [ \%ind0, \%ind1, ... ]   # one hashref per sample/individual
#   each %ind = ( marker_name => [allele, allele, ...] )   # ploidy alleles per marker
#   Alleles equal to 'N' are treated as MISSING and excluded from per-site sample size.
use strict; 
use warnings; 

# {allele=>count} for one marker across all individuals, excluding 'N'.
sub _site_counts {
	my ($inds, $m) = @_; 
	my %c; 
	for my $ind (@$inds) {
		for my $a ( @{ $ind->{$m} // [] } ) {
			$a eq 'N' and next; 
			$c{$a}++; 
		}
	}
	return \%c; 
}

sub marker_names { return ( keys %{ $_[0][0] } ); }

# number of segregating (polymorphic, >1 allele) sites; N excluded.
sub segregating_sites_count {
	my ($inds) = @_; 
	my $S = 0; 
	for my $m ( marker_names($inds) ) {
		scalar( keys %{ _site_counts($inds, $m) } ) > 1 and $S++; 
	}
	return $S; 
}

# nucleotide diversity pi = sum over sites of avg pairwise differences (per-site n, N excluded).
sub pi {
	my ($inds) = @_; 
	my $pi = 0; 
	for my $m ( marker_names($inds) ) {
		my $c = _site_counts($inds, $m); 
		my ($n, $sum2) = (0, 0); 
		for my $cnt ( values %$c ) { $n += $cnt; $sum2 += $cnt*$cnt; }
		$n < 2 and next; 
		$pi += ($n*$n - $sum2) / ( $n * ($n-1) );   # = n/(n-1) * (1 - sum p_i^2)
	}
	return $pi; 
}

# Watterson theta = sum over segregating sites of 1/a1(n_site); N excluded (matches _selfTheta intent).
sub theta {
	my ($inds) = @_; 
	my $theta = 0; 
	for my $m ( marker_names($inds) ) {
		my $c = _site_counts($inds, $m); 
		my $n = 0; $n += $_ for values %$c; 
		$n < 2 and next; 
		scalar( keys %$c ) > 1 or next; 
		my $a1 = 0; $a1 += 1/$_ for (1 .. $n-1); 
		$a1 > 0 and $theta += 1/$a1; 
	}
	return $theta; 
}

# Tajima's D over a window. Uses a single nominal sample size $n (default: the
# largest per-site sample size in the window = full sample when no missing data).
sub tajima_D {
	my ($inds, $n) = @_; 
	my $S = segregating_sites_count($inds); 
	$S > 0 or return undef; 
	unless ( defined $n ) {
		$n = 0; 
		for my $m ( marker_names($inds) ) {
			my $t = 0; $t += $_ for values %{ _site_counts($inds, $m) }; 
			$t > $n and $n = $t; 
		}
	}
	$n > 1 or return undef; 
	my ($a1, $a2) = (0, 0); 
	for my $i (1 .. $n-1) { $a1 += 1/$i; $a2 += 1/($i*$i); }
	my $b1 = ($n+1) / ( 3*($n-1) ); 
	my $b2 = 2*($n*$n + $n + 3) / ( 9*$n*($n-1) ); 
	my $c1 = $b1 - 1/$a1; 
	my $c2 = $b2 - ($n+2)/($a1*$n) + $a2/($a1*$a1); 
	my $e1 = $c1 / $a1; 
	my $e2 = $c2 / ($a1*$a1 + $a2); 
	my $Vd = $e1*$S + $e2*$S*($S-1); 
	$Vd > 0 or return undef; 
	my $pi_val = pi($inds); 
	my $theta_w = $S / $a1; 
	return ( $pi_val - $theta_w ) / sqrt($Vd); 
}

1; 
