package SeqAlnSunhh; 

#Author: Honghe Sun, hs738@cornell.edu, 2014-07-15
# 2014-07-15 In function olap_e2e_A2B(), "N" is always treated as "mismatch"; 

use strict; 
use warnings; 
use LogInforSunhh; 
use Exporter qw(import); 

our @EXPORT = qw(olap_e2e_A2B); 
our @EXPORT_OK; 

# define the scores for match, mismatch and indel:
my %para; 
$para{match_score} = 2; 
$para{mismatch_score} = -1; 
$para{gapopen_score} = 0; 
$para{gapextend_score} = -2; 
$para{opt_tail} = 0; 
$para{opt_head} = 0; 

## Do not distinguish "N" or other special words. 
## Input  : (sequenceA, sequenceB [, {%para}])
## Output : { qw/(start|end)(A|B) aln_len aln_ident count seqA_aln seqB_aln seqC_aln/ => values }
##          startA, endA : start and end positions of seqA aligned; 
##          startB, endB : start and end positions of seqB aligned; 
##          aln_len      : length of aligned region (including gap)
##          aln_ident    : (length of match bases) / aln_len
##          seqA_aln     : alignment in seqA
##          seqB_aln     : alignemtn in seqB
##          seqC_aln     : Consensus signs between seqA_aln and seqB_aln. "." - same, "D" - diff, "-" - gap. 
sub olap_e2e_A2B {
	my ($seqA, $seqB, $ref_para) = @_; 
	my %loc_para = %para; 
	if (defined $ref_para and ref($ref_para) eq 'HASH') {
		for my $tk ( keys %$ref_para ) {
			$loc_para{$tk} = $ref_para->{$tk}; 
		}
	}
	
	# return back values; 
	my ($psA, $peA, $psB, $peB); # (startA, endA, startB, endB) of overlap region. 
	my $seqA_aln = ''; # Aligned sequence from seqA 
	my $seqB_aln = ''; # Aligned sequence from seqB 
	my ($aln_len, $aln_ident) = (0, 0); 

	# Caculating. 
	$seqA = uc($seqA); 
	$seqB = uc($seqB); 
	my $lenA = length($seqA); 
	my $lenB = length($seqB); 
	my @baseA = ( $seqA =~ m/(.)/g ); 
	my @baseB = ( $seqB =~ m/(.)/g ); 
	my @smat;               # array matrix storing the dynamic programming score matrix; 
	my @traceback;          # array storing the traceback directions. 
	
	# Define two-dimensional matrices "@smat" and "@traceback" with ($lenA + 1) rows, and ($lenB + 1) columns, and initialize it with 0s. 
	for (my $i=0; $i<=$lenA; $i++) {
		for (my $j=0; $j<=$lenB; $j++) {
			$smat[$i][$j]      = 0; 
			$traceback[$i][$j] = 0; 
		}
	}
	
	# Setup boundary of traceback array, used to stop searching. 
	# Put "."s in the first row of matrix @traceback; 
	for (my $j=0; $j<=$lenB; $j++) {
		$traceback[0][$j] = "."; 
	}
	# Put "."s in the first column of matrix @traceback; 
	for (my $i=0; $i<=$lenA; $i++) {
		$traceback[$i][0] = "."; 
	}
	
	# Dynamic programming recursion
	for (my $i=1; $i<=$lenA; $i++) {
		# base ($baseA[$i-1]) position $i in $seqA (seqA goes down the first column) 
		for (my $j=1; $j<=$lenB; $j++) {
			# base ($baseB[$j-1]) position $j in $seqB (seqB goes across the first row)
			my ($score , $diag, $up, $left, $max) ; # scores used to construct matrix; 
			
			# Find the value to put into $smat[$i][$j]
			# (1) the first possibility is to take the diagonal element plus a match/mismatch score; 
			if ( $baseA[$i-1] eq $baseB[$j-1] and $baseA[$i-1] ne "N" ) { $score = $loc_para{match_score};    } 
			else                                { $score = $loc_para{mismatch_score}; }
			$diag = $smat[$i-1][$j-1]+$score; 
			# (2) the second possibility is to take the element above plus a gap score; 
			if ( $traceback[$i-1][$j  ] eq '|' ) { $up = $smat[$i-1][$j  ] + $loc_para{gapextend_score};                        } 
			else                                 { $up = $smat[$i-1][$j  ] + $loc_para{gapextend_score} + $loc_para{gapopen_score}; }
			# (3) the third possibility is to take the element on the left plus a gap score; 
			if ( $traceback[$i  ][$j-1] eq '-' ) { $left = $smat[$i  ][$j-1] + $loc_para{gapextend_score};                        }
			else                                 { $left = $smat[$i  ][$j-1] + $loc_para{gapextend_score} + $loc_para{gapopen_score}; }
			# record which of the three possibilities was highest, into @traceback: 
			if    ( $diag >  $up   && $diag >  $left ) { $traceback[$i][$j] = '>'; $max = $diag; }
			elsif ( $up   >  $diag && $up   >  $left ) { $traceback[$i][$j] = '|'; $max = $up;   }
			elsif ( $left >  $diag && $left >  $up   ) { $traceback[$i][$j] = '-'; $max = $left; }
			elsif ( $left == $diag && $up   == $diag ) { $traceback[$i][$j] = '*'; $max = $diag; }
			elsif ( $up   == $left && $up   >  $diag && $left >  $diag ) { $traceback[$i][$j] = 'L'; $max = $up;   }
			elsif ( $up   == $diag && $up   >  $left && $diag >  $left ) { $traceback[$i][$j] = 'V'; $max = $up;   }
			elsif ( $diag == $left && $diag >  $up   && $left >  $up   ) { $traceback[$i][$j] = 'Z'; $max = $diag; }
			else { &stopErr( "Why we come here!\n" );  }
			# record the highest score in @smat; 
			$smat[$i][$j] = $max; 
		}
	}
	
	# The best overlap is the given by the largest number in the bottom-most row of the matrix @smat; 
	# Because we are doing A-2-B overlap alignment; 
	my $best_x = -1; # Position aligned in seqA 
	my $best_y = -1; # Position aligned in seqB 
	my $best_score = -1e15; 
	
	# Here we only check the end of seqA, and we want to get the left-most alignment for A-2-B end-to-end overlap 
	# So we check only the bottom-most row from left to right. 
	if ($loc_para{opt_tail} > 0) {
		for (my $i=$lenA; $i>=0 && $i >= $lenA-$loc_para{opt_tail} ; $i--) {
			for (my $j=0; $j<=$lenB; $j++) {
				if ( $smat[$i][$j] > $best_score ) { $best_x = $i; $best_y = $j; $best_score = $smat[$i][$j]; } 
			}
		}
	}else{
		for (my $j=0; $j<=$lenB; $j++) { 
			# If we want right-most alignment on seqB, we should use "for (my $j=$lenB; $j>=0; $j--)"
			if ( $smat[$lenA][$j] > $best_score ) { $best_x = $lenA; $best_y = $j; $best_score = $smat[$lenA][$j]; } 
		}
	}
	
	# Record the boundary of best overlap region. 
	$peA = $best_x; 
	$peB = $best_y; 
	
	## record the traceback as '#': 
	# construct the best alignment. 
	my $tracebackvalue = $traceback[$best_x][$best_y]; 
	my $prev_score = $smat[$best_x][$best_y]; 
	my %count = qw(same 0 diff 0 gapA 0 gapB 0); 
	my $seqC_aln = ''; 
	while ( $tracebackvalue ne '.' ) {
		$traceback[$best_x][$best_y] = "#$traceback[$best_x][$best_y]"; 
		my ($prev_x, $prev_y) = ($best_x, $best_y); 
		my ($seqA_letter, $seqB_letter, $seqC_letter, $cKey); 
		if ( $tracebackvalue eq '>' or $tracebackvalue eq '*' or $tracebackvalue eq 'V' or $tracebackvalue eq 'Z' ) {
			$seqA_letter = $baseA[$best_x-1]; 
			$seqB_letter = $baseB[$best_y-1]; 
			if ( $seqA_letter eq $seqB_letter ) {
				$cKey = 'same'; 
				$seqC_letter = '.'; 
			} else {
				$cKey = 'diff'; 
				$seqC_letter = 'D'; 
			}
			$best_x --; 
			$best_y --; 
		} elsif ( $tracebackvalue eq '-' or $tracebackvalue eq 'L' ) {
			$seqA_letter = '-'; 
			$seqB_letter = $baseB[$best_y-1]; 
			$cKey = 'gapA'; 
			$seqC_letter = '-'; 
			$best_y --; 
		} elsif ( $tracebackvalue eq '|' ) {
			$seqA_letter = $baseA[$best_x-1]; 
			$seqB_letter = '-'; 
			$cKey = 'gapB'; 
			$seqC_letter = '-'; 
			$best_x --; 
		} else {
			&tsmsg( "[Err] Unknown tracebackvalue [$tracebackvalue]\n" ); 
			exit 1; 
		}
		$seqA_aln = $seqA_letter . $seqA_aln; 
		$seqB_aln = $seqB_letter . $seqB_aln; 
		$seqC_aln = $seqC_letter . $seqC_aln; 
		$count{$cKey} ++; 
		$aln_len ++; 
		if ( $loc_para{opt_head} > 0 ) {
			if ( $best_y <= $loc_para{opt_head} and $prev_score < $smat[$best_x][$best_y] ) {
				# ($best_x, $best_y) = ($prev_x, $prev_y); 
				last; 
			}
		}
		# Find the new traceback value; 
		$tracebackvalue = $traceback[$best_x][$best_y]; 
		$prev_score = $smat[$best_x][$best_y]; 
	}
	$psA = $best_x+1; 
	$psB = $best_y+1; 
	
	## Prepare to output score matrix (@smat) and traceback matrix (@traceback)
	## Put $seqB in the first row of matrix @smat; 
	for (my $j=1; $j<=$lenB; $j++) {
		$smat[0][$j] = $baseB[$j-1]; 
	}
	## Put $seqA in the first column of matrix @smat; 
	for (my $i=1; $i<=$lenA; $i++) {
		$smat[$i][0] = $baseA[$i-1]; 
	}
#	## Print out the result. 
#	print "Overlap matrix, with traceback shown as #:\n"; 
#	my $smatvalue; 
#	$tracebackvalue = undef(); 
#	for (my $i=0; $i<=$lenA; $i++) {
#		for (my $j=0; $j<=$lenB; $j++) {
#			$smatvalue = $smat[$i][$j]; 
#			$tracebackvalue = $traceback[$i][$j]; 
#			print "$smatvalue $tracebackvalue\t"; 
#		}
#		print "\n"; 
#	}
#	print "\n"; 
	
#	# Print out the overlap alignment; 
#	print "Overlap alignment (score = $best_score):\n"; 
#	print "A: $seqA_aln\n"; 
#	print "B: $seqB_aln\n"; 
#	print "C: $seqC_aln\n"; 
	
	$aln_ident = ( $aln_len == 0 ) ? -1 : ($count{same})/$aln_len ; 
	my %back; 
	$back{startA} = $psA; $back{endA} = $peA; 
	$back{startB} = $psB; $back{endB} = $peB; 
	$back{aln_len} = $aln_len; $back{aln_ident} = $aln_ident; 
	$back{aln_score} = $best_score; 
	$back{count} = \%count; 
	$back{seqA_aln} = $seqA_aln; $back{seqB_aln} = $seqB_aln; 
	$back{seqC_aln} = $seqC_aln; 
	return \%back; 
}#End sub find_overlap_alignment

1; # It is important to include this line. 

