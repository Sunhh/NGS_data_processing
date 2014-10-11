#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 dgt.tab\n"; 

# eleID   eleS    eleE    Str     seqID   LTR1_S  LTR1_E  LTR2_S  LTR2_E  Inner_S Inner_E PBS_S   PBS_E   PPT_S   PPT_E   scfID
# RR1     98334   105001  ?       seq10   98339   99198   104137  104996  99199   104136  -1      -1      -1      -1      S400016_pilon
# RR2     785204  788754  +       seq10   785209  785369  788589  788749  785370  788588  785373  785384  -1      -1      S400016_pilon

my $max_dist_to_ltr = 20; 
my $minR_inInclRegion = 0.5; 

my %rec; 
my @id_list; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[0] eq 'eleID' ) {
		print STDOUT "$_\n"; 
		next; 
	}
	my ($inner_s, $inner_e, $pbs_s, $pbs_e, $ppt_s, $ppt_e) = @ta[9,10, 11,12, 13,14]; 
	my ($eleID, $seqID, $scfID) = @ta[0,4,15]; 
	my $pbs_olap = &olapLen( $inner_s, $inner_e, $pbs_s, $pbs_e ); 
	my $ppt_olap = &olapLen( $inner_s, $inner_e, $ppt_s, $ppt_e ); 
	my $pbs_dist = &minDist( $inner_s, $inner_e, $pbs_s, $pbs_e ); 
	my $ppt_dist = &minDist( $inner_s, $inner_e, $ppt_s, $ppt_e ); 
	my $is_good = 0; 
	&good_lt( $pbs_dist, $max_dist_to_ltr ) and &good_gt( $pbs_olap, $minR_inInclRegion * ($pbs_e-$pbs_s+1) ) and $is_good = 1; 
	&good_lt( $ppt_dist, $max_dist_to_ltr ) and &good_gt( $ppt_olap, $minR_inInclRegion * ($ppt_e-$ppt_s+1) ) and $is_good = 1; 
	$is_good == 1 or next; 
	print STDOUT "$_\n"; 
	#my $tk1 = "${eleID}_${seqID}_$ta[1]_$ta[2]_ELE_$scfID"; # Element region. (with TSD)
	#my $tk2 = "${eleID}_${seqID}_$ta[5]_$ta[8]_LTR_$scfID"; # LTR region. (without TSD)
	#my $tk3 = "${eleID}_${seqID}_$ta[9]_$ta[10]_INN_$scfID"; # Internal region. (without ltr region)
	#print STDOUT "$tk1\n$tk2\n$tk3\n"; 
}

sub good_lt {
	my ($v1, $v2) = @_; 
	$v1 != -1 and $v1 <= $v2 and return 1; 
	return 0; 
}
sub good_gt {
	my ($v1, $v2) = @_; 
	$v1 != -1 and $v1 >= $v2 and return 1; 
	return 0; 
}

sub olapLen {
	my ($s1,$e1, $s2, $e2) = @_; 
	( $s1 < 0 or $s2 < 0 ) and return (-1); 
	my $s3 = ($s1 > $s2) ? $s1 : $s2 ; 
	my $e3 = ($e1 < $e2) ? $e1 : $e2 ; 
	return ($e3-$s3+1) ; 
}

sub minDist {
	my ($s1,$e1, $s2, $e2) = @_; 
	( $s1 < 0 or $s2 < 0 ) and return (-1); 
	$e2 < $s1 and do { &tsmsg("[Err] Skip minDist [$e2 < $s1 : ($s1,$e1, $s2, $e2)].\n"); return (-1); }; 
	$s2 > $e1 and do { &tsmsg("[Err] Skip minDist [$s2 > $e1 : ($s1,$e1, $s2, $e2)].\n"); return (-1); }; 
	my $dist1 = $s2-$s1; $dist1 < 0 and $dist1 = 0; 
	my $dist2 = $e1-$e2; $dist2 < 0 and $dist2 = 0; 
	my $dist = ($dist1 > $dist2) ? $dist2  : $dist1 ; 
	return $dist ; 
}



