#!/usr/bin/perl
use strict;
use warnings; 
use ReadInAlnSunhh; 
use Getopt::Long; 

my %opts; 
GetOptions(\%opts, 
	"help!", 
	"minClustLen:i", 
	"maxEdgeDist:i", 
	"minRatioShrt:f", 
	"minRatioIncl:f", 
	"out:s", 
); 

sub usage {
	print STDOUT <<HELP; 
perl $0 in_syn.maf
-help    
-minClustLen      [500] Minimum length of clustered syn_blocks; 
-maxEdgeDist      [10000] Maximum length of distance between terminal block and scaffold edges. 
-minRatioShrt     [0.9] Minimum ratio of short scaffolds (<= maxEdgeDist) required for a good alignment. 
-minRatioIncl     [0.9] Minimum ratio of included scaffolds in one block. 
HELP
	exit 1; 
}
-t and !@ARGV and &usage(); 
$opts{help} and &usage(); 

defined $opts{minClustLen} or $opts{minClustLen} = 500; 
defined $opts{maxEdgeDist} or $opts{maxEdgeDist} = 10000; 
defined $opts{minRatioShrt} or $opts{minRatioShrt} = 0.9; # Maybe 0.6 will be better, but I am not sure about it. 
defined $opts{minRatioIncl} or $opts{minRatioIncl} = 0.9; 

my $oFh = \*STDOUT; 
if (defined $opts{out}) {
	my $tfh; 
	open $tfh,'>',"$opts{out}" or die "1: $!\n"; 
	$oFh = $tfh; 
}
my @FH; 
!(-t) and push(@FH, \*STDIN); 
for (@ARGV) {
	my $fh; 
	open $fh, '<', "$_" or die "2: $!\n"; 
	push(@FH, $fh); 
}

my %all_links; #
my %rep_incl; 
my %seqLength; 
my (%r2q, %q2r); 
# key1 = Ref_scaf_id 
#  key2 = {p5, p3, incl, cont, same} # { "5p extend to", "3p extend to", "wholely included by", "contain whole sequence of", "same to" }
#   Value of "p5"   : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR] }, ...
#   Value of "p3"   : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR] }, ...
#   Value of "incl" : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR] }, ...
#   Value of "cont" : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR] }, ...
#   Value of "same" : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR] }, ...


&tsmsg("[Rec] Reading in MAF files.\n"); 
for my $fh (@FH) {
	while ( my %rec1 = %{readMAF($fh)} ) {
		# Check if it is a pairwise-alignment . 
		@{$rec1{o}} >= 2 or next; 
		my @sLines; 
		for my $tl (@{$rec1{o}}) {
			$tl =~ m/^s\s/ and push(@sLines, $tl); 
		}
		@sLines >= 2 or next; 
		
		# Get Ref/Qry alignment information. 
		my %s1 = %{ splitMafSline($sLines[0], 1) }; # Reference alignment 
		my %s2 = %{ splitMafSline($sLines[1], 1) }; # Query     alignment 
		my $tmp_str = 'f'; $s1{seqStrand} ne $s2{seqStrand} and $tmp_str = 'r'; 
		defined $seqLength{ $s1{seqId} } or $seqLength{ $s1{seqId} } = $s1{seqLen}; 
		defined $seqLength{ $s2{seqId} } or $seqLength{ $s2{seqId} } = $s2{seqLen}; 
		
		# Check if this alignment fits the requirements. 
		if ( $s1{seqLen} >= $opts{minClustLen} and $s2{seqLen} >= $opts{minClustLen} ) {
			( $s1{blkSize} >= $opts{minClustLen} and $s2{blkSize} >= $opts{minClustLen} ) or next; # Both aligned blocks should be larger than minimum clustered length 
		} else {
			if ( $s1{seqLen} < $opts{minClustLen} ) { $s1{blkSize} >= $s1{seqLen} * $opts{minRatioShrt} or next;  } 
			if ( $s2{seqLen} < $opts{minClustLen} ) { $s2{blkSize} >= $s2{seqLen} * $opts{minRatioShrt} or next;  } 
		}
		if ( $s1{seqLen} > $s2{seqLen} ) {
			$s2{seqLen} <= $opts{maxEdgeDist} and do { $s2{blkSize} >= $opts{minRatioShrt} * $s2{seqLen} or next; }; 
		} else {
			$s1{seqLen} <= $opts{maxEdgeDist} and do { $s1{blkSize} >= $opts{minRatioShrt} * $s1{seqLen} or next; }; 
		}
		
		# Check if this alignment is in the end of a scaffold. 
		## Store links that maybe useful. 
		## Discard the seqLen of scaffolds now. 
		if ( $s1{normS} <= $opts{maxEdgeDist} ) {
			if ( $s1{normE} >= $s1{seqLen} - $opts{maxEdgeDist} + 1 ) {
				# [incl] s1 is included by Qry (s2); Maybe s1 is same to s2, but I don't care here. 
				if ( defined $all_links{$s1{seqId}}{incl}{$s2{seqId}} ) { 
					# defined Repeat included information for s1 - s2 pair. This should be impossible. 
					&tsmsg("[Wrn] Repeat linkage [included] between [$s1{seqId}] and [$s2{seqId}]\n");
					## Now there are multiple [included] blocks, so these blocks should be separated into [p5/p3] types. 
					## We should remove these [included] bocks before searching links. 
					$rep_incl{$s1{seqId}}{$s2{seqId}} = 1; 
					defined $all_links{$s1{seqId}}{p5}{$s2{seqId}} or $all_links{$s1{seqId}}{p5}{$s2{seqId}} = [ @{ $all_links{$s1{seqId}}{incl}{$s2{seqId}} } ]; 
					$s1{normS} < $all_links{$s1{seqId}}{p5}{$s2{seqId}}[5] and $all_links{$s1{seqId}}{p5}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
					defined $all_links{$s1{seqId}}{p3}{$s2{seqId}} or $all_links{$s1{seqId}}{p3}{$s2{seqId}} = [ @{ $all_links{$s1{seqId}}{incl}{$s2{seqId}} } ]; 
					$s1{normE} > $all_links{$s1{seqId}}{p3}{$s2{seqId}}[6] and $all_links{$s1{seqId}}{p3}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ];

					$s1{blkSize} > $all_links{$s1{seqId}}{incl}{$s2{seqId}}[7] and $all_links{$s1{seqId}}{incl}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
				} else {
					$all_links{$s1{seqId}}{incl}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
				}
				$q2r{$s2{seqId}}{$s1{seqId}} = 1; 
			} else {
				# [p5] s1 should be on the 3p side of s2 , so s1 is extended to 5p side according to s2. 
				if ( defined $all_links{$s1{seqId}}{p5}{$s2{seqId}} ) {
					if ( $s1{normS} < $all_links{$s1{seqId}}{p5}{$s2{seqId}}[5] ) {
						&tsmsg("[Wrn] More to-terminal linkage [p5] between [$s1{seqId}] and [$s2{seqId}]\n"); 
						$all_links{$s1{seqId}}{p5}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
					} else {
						; 
					}
				} else {
					$all_links{$s1{seqId}}{p5}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
				}
				$q2r{$s2{seqId}}{$s1{seqId}} = 1; 
			}# End if ( $s1{normE} >= $s1{seqLen} - $opts{maxEdgeDist} + 1 )
		} elsif ( $s1{normE} >= $s1{seqLen} - $opts{maxEdgeDist} + 1 ) { 
			# [p3] s1 should be extended to 3p side according to s2; 
			$q2r{$s2{seqId}}{$s1{seqId}} = 1; 
			if ( defined $all_links{$s1{seqId}}{p3}{$s2{seqId}} ) {
				if ( $s1{normE} > $all_links{$s1{seqId}}{p3}{$s2{seqId}}[6] ) {
					&tsmsg("[Wrn] More to-terminal linkage [p3] between [$s1{seqId}] and [$s2{seqId}]\n"); 
					$all_links{$s1{seqId}}{p3}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
				} else {
					; 
				}
			} else {
				$all_links{$s1{seqId}}{p3}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
			}# End if ( defined $all_links{$s1}{p3}{$s2{seqId}} )
		} elsif ( $s2{normS} <= $opts{maxEdgeDist} and $s2{normE} >= $s2{seqLen} - $opts{maxEdgeDist} + 1 ) { 
			# [cont] s1 should contain the whole sequence of s2. However, this information is not useful for linking scaffolds from Ref. 
			$q2r{$s2{seqId}}{$s1{seqId}} = 1; 
			if ( defined $all_links{$s1{seqId}}{cont}{$s2{seqId}} ) {
				&tsmsg("[Wrn] Repeat containing linkage [cont] between [$s1{seqId}] and [$s2{seqId}]\n"); 
				$s2{blkSize} > $all_links{$s1{seqId}}{cont}{$s2{seqId}}[4] and $all_links{$s1{seqId}}{cont}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
			} else {
				$all_links{$s1{seqId}}{cont}{$s2{seqId}} = [ $s2{seqLen}, $tmp_str, @s2{qw/normS normE blkSize/}, @s1{qw/normS normE blkSize/}, $s2{seqId}, $s1{seqId} ]; 
			}# End if ( defined $all_links{$s1}{cont}{$s2{seqId}} )
		} else {
			# An internal linkage that is not used. 
			; 
		}
	}
}#End for my $fh (@FH) 

#############################################################
# Checking 
{
# Remove block links in %rep_incl
&tsmsg("[Rec] Deleting repeated [included] links which have been added to [p5/p3] links.\n"); 
my $cnt = 0; 
for my $incl1 ( keys %rep_incl ) {
	for my $incl2 ( keys %{$rep_incl{$incl1}} ) {
		delete $all_links{$incl1}{incl}{$incl2} ; 
		$cnt ++; 
	}
}
&tsmsg("[Rec] [$cnt] links deleted.\n"); 
# Remove [incl] block links lower than $opts{minRatioIncl}
my %shrt_incl; 
for my $incl1 (keys %all_links) {
	defined $all_links{$incl1}{incl} or next; 
	for my $incl2 (keys %{$all_links{$incl1}{incl}}) {
		$all_links{$incl1}{incl}{$incl2}[7] >= $opts{minRatioIncl} * $seqLength{$incl1} or $shrt_incl{$incl1}{$incl2} = 1; 
	}
}
$cnt = 0; 
for my $incl1 (keys %shrt_incl) {
	for my $incl2 (keys %{$shrt_incl{$incl1}}) {
		delete $all_links{$incl1}{incl}{$incl2}; 
		$cnt ++; 
	}
}
&tsmsg("[Rec] [$cnt] links deleted.\n"); 
}


#############################################################
# Tracing links
&tsmsg("[Rec] Tracing good links.\n"); 
my %good_links; 
for my $refId ( keys %all_links ) {
	if ( defined $all_links{$refId}{p5} ) {
		# Ref_link_type -- "p5"
		my @qIDs = keys %{$all_links{$refId}{p5}}; 
		if ( scalar( @qIDs ) == 1 ) {
			# $refId may be extended to 5p side. 
			my $qryId = $qIDs[0]; 
			my $nearest_link; # Same to element of all_links; 
			for my $r2Id ( keys %{$q2r{$qryId}} ) {
				$r2Id eq $refId and next; 
				for my $r2_linkType ( qw/p5 p3 incl/ ) {
					( defined $all_links{$r2Id}{$r2_linkType} and defined $all_links{$r2Id}{$r2_linkType}{$qryId} ) or next; 
					# ( keys %{ $all_links{$r2Id}{$r2_linkType} } ) == 1 or next; 
					$nearest_link = &closer_r2q_p5(
					  $all_links{$refId}{p5}{$qryId}, 
					  $all_links{$r2Id}{$r2_linkType}{$qryId}, 
					  $r2_linkType, 
					  $nearest_link
					); 
				}
			}#End for my $r2Id ( keys %q2r{$qryId} )
			&save_good_link($all_links{$refId}{p5}{$qryId}, $nearest_link, \%good_links); 
		}
	} elsif ( defined $all_links{$refId}{p3} ) {
		# Ref_link_type -- "p3"
		my @qIDs = keys %{$all_links{$refId}{p3}}; 
		if ( scalar( @qIDs ) == 1 ) {
			# $refId may be extended to 3p side. 
			my $qryId = $qIDs[0]; 
			my $nearest_link; 
			for my $r2Id ( keys %{$q2r{$qryId}} ) {
				$r2Id eq $refId and next; 
				for my $r2_linkType ( qw/p5 p3 incl/ ) {
					( defined $all_links{$r2Id}{$r2_linkType} and defined $all_links{$r2Id}{$r2_linkType}{$qryId} ) or next; 
					# ( keys %{ $all_links{$r2Id}{$r2_linkType} } ) == 1 or next; 
					$nearest_link = &closer_r2q_p3(
					  $all_links{$refId}{p3}{$qryId}, 
					  $all_links{$r2Id}{$r2_linkType}{$qryId}, 
					  $r2_linkType, 
					  $nearest_link
					); 
				}
			}#End for my $r2Id ( keys %q2r{$qryId} )
			&save_good_link($all_links{$refId}{p3}{$qryId}, $nearest_link, \%good_links); 
		}
	} elsif ( defined $all_links{$refId}{incl} ) { 
		# Ref_link_type -- "incl"
		my @qIDs = keys %{$all_links{$refId}{incl}}; 
		if ( scalar( @qIDs ) == 1 ) {
			# $refId may be extended to either 5p/3p side. 
			## For 5p/3p side
			my $qryId = $qIDs[0]; 
			my $nearest_link_p5; # Same to element of all_links; 
			my $nearest_link_p3; 
			for my $r2Id ( keys %{$q2r{$qryId}} ) {
				$r2Id eq $refId and next; 
				for my $r2_linkType ( qw/p5 p3 incl/ ) {
					( defined $all_links{$r2Id}{$r2_linkType} and defined $all_links{$r2Id}{$r2_linkType}{$qryId} ) or next; 
					# ( keys %{ $all_links{$r2Id}{$r2_linkType} } ) == 1 or next; 
					$nearest_link_p5 = &closer_r2q_p5(
					  $all_links{$refId}{incl}{$qryId}, 
					  $all_links{$r2Id}{$r2_linkType}{$qryId}, 
					  $r2_linkType, 
					  $nearest_link_p5
					); 
					$nearest_link_p3 = &closer_r2q_p3(
					  $all_links{$refId}{incl}{$qryId}, 
					  $all_links{$r2Id}{$r2_linkType}{$qryId}, 
					  $r2_linkType, 
					  $nearest_link_p3
					); 
				}
			}#End for my $r2Id ( keys %q2r{$qryId} )
			&save_good_link($all_links{$refId}{incl}{$qryId}, $nearest_link_p5, \%good_links); 
			&save_good_link($all_links{$refId}{incl}{$qryId}, $nearest_link_p3, \%good_links); 
		}
	} else {
		; 
	}# End 'else' from if ( defined $all_links{$refId}{p5} ) 
}

&tsmsg("[Rec] Output good links.\n"); 
print {$oFh} join("\t", qw/Ref1ID Ref2ID Ref1Len Ref2Len Ref1_blk_normS Ref1_blk_normE Ref1_blk_fr Ref2_blk_normS Ref2_blk_normE Ref2_blk_fr Gap_bwt_Blocks QryID QryLen QryBlk1 QryBlk2/)."\n"; 
for my $refId (sort keys %good_links) {
	for my $qryId (sort keys %{$good_links{$refId}}) {
		my @ta1 = @{$good_links{$refId}{$qryId}}; 
		print {$oFh} join("\t", $refId, $qryId, $seqLength{$refId}, $seqLength{$qryId}, @ta1)."\n"; 
	}
}
&tsmsg("[Rec] All done.\n"); 

# my %all_links; #
# key1 = Ref_scaf_id 
#  key2 = {p5, p3, incl, cont, same} # { "5p extend to", "3p extend to", "wholely included by", "contain whole sequence of", "same to" }
#   Value of "p5"   : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR, Qry_scaf_id, Ref_scaf_id] }, ...
#   Value of "p3"   : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR, Qry_scaf_id, Ref_scaf_id], }, ...
#   Value of "incl" : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR, Qry_scaf_id, Ref_scaf_id], }, ...
#   Value of "cont" : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR, Qry_scaf_id, Ref_scaf_id], }, ...
#   Value of "same" : {"Qry_scaf_id" => [Len_Qry, "+/-", Start_Position_Qry, End_Position_Qry, blkSizeQ, Start_Position_Ref, End_Position_Ref, blkSizeR, Qry_scaf_id, Ref_scaf_id], }, ...

## Output : Change %good_links : {refId}{qryId} => [ Ref_blk_normS, Ref_blk_normE, Ref_blk_fr, Qry_blk_normS, Qry_blk_normE, Qry_blk_fr ]
sub save_good_link {
	my ( $r0_link, $r1_link, $save_link) = @_; 
	( defined $r1_link and ref($r1_link) eq 'ARRAY' and @$r1_link > 0) or return 1; 
	my @v0 = @$r0_link; 
	my @v1 = @$r1_link; 
	$v0[8] eq $v1[8] or die "3: $!\nv0:@v0\nv1:@v1\n"; 
	if ( defined $save_link->{$v0[9]}{$v1[9]} ) {
		&tsmsg("[Msg] Multiple good links between [$v0[9]] and [$v1[9]] :\nV0:@v0\nV1:@v1\n"); 
	} else {
		if ( $v0[3] <= $v1[2] ) {
			# s0 -> s1 
			my $gapLen = $v1[2]-$v0[3]-1; 
			if ( $v0[1] eq 'f' ) {
				if ( $v1[1] eq 'f' ) {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'f', @v1[5,6], 'f', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				} else {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'f', @v1[5,6], 'r', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				}
			} else {
				if ( $v1[1] eq 'f' ) {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'r', @v1[5,6], 'f', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				} else {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'r', @v1[5,6], 'r', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				}
			}
		} elsif ( $v0[2] >= $v1[3] ) {
			# s1 -> s0
			my $gapLen = $v0[2]-$v1[3]-1; 
			if ( $v0[1] eq 'f' ) {
				if ( $v1[1] eq 'f' ) {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'r', @v1[5,6], 'r', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				} else {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'r', @v1[5,6], 'f', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				}
			} else {
				if ( $v1[1] eq 'f' ) {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'f', @v1[5,6], 'r', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				} else {
					$save_link->{$v0[9]}{$v1[9]} = [ @v0[5,6], 'f', @v1[5,6], 'f', $gapLen, $v0[8], $v0[0], "$v0[1]:$v0[2]-$v0[3]", "$v1[1]:$v1[2]-$v1[3]" ]; 
				}
			}
		} else {
			&tsmsg("[Wrn] bad $v0[3] VS. $v1[2]\n"); 
		}
	}
	return 0; 
}


sub closer_r2q_p5 {
	my ($r0_link, $r1_link, $r1_link_type, $prev_link) = @_; 
	my $back_link = {}; 
	my ($r0_str, $r0_QS, $r0_QE) = (@{$r0_link}[1,2,3]); 
	my ($r1_str, $r1_QS, $r1_QE) = (@{$r1_link}[1,2,3]); 
	my ($prev_str, $prev_QS, $prev_QE); 
	if (defined $prev_link and ref($prev_link) eq 'ARRAY' and @$prev_link > 0) {
		($prev_str, $prev_QS, $prev_QE) = (@{$prev_link}[1,2,3]); 
		$back_link = $prev_link; 
	} else {
		($prev_str, $prev_QS, $prev_QE) = ('', '', ''); 
	}
	if ( 
	  ( $r1_link_type eq 'p5' and $r1_str ne $r0_str ) 
	  or
	  ( $r1_link_type eq 'p3' and $r1_str eq $r0_str ) 
	  or
	  ( $r1_link_type eq 'incl'                      ) 
	) {
		if ( $r0_str eq 'f' ) {
			&left_closer($r0_QS, $r0_QE, $r1_QS, $r1_QE, $prev_QS, $prev_QE) == 1 and $back_link = $r1_link; 
		} elsif ( $r0_str eq 'r' ) {
			&right_closer($r0_QS, $r0_QE, $r1_QS, $r1_QE, $prev_QS, $prev_QE) == 1 and $back_link = $r1_link; 
		} else {
			die "4: $!\n|$r0_str|\n"; 
		}
	} else {
		return $back_link; 
	}
	
	return $back_link; 
}# End sub closer_r2q_p5()

sub closer_r2q_p3 {
	my ($r0_link, $r1_link, $r1_link_type, $prev_link) = @_; 
	my $back_link = {}; 
	my ($r0_str, $r0_QS, $r0_QE) = (@{$r0_link}[1,2,3]); 
	my ($r1_str, $r1_QS, $r1_QE) = (@{$r1_link}[1,2,3]); 
	my ($prev_str, $prev_QS, $prev_QE); 
	if (defined $prev_link and ref($prev_link) eq 'ARRAY' and @$prev_link > 0) {
		($prev_str, $prev_QS, $prev_QE) = (@{$prev_link}[1,2,3]); 
		$back_link = $prev_link; 
	} else {
		($prev_str, $prev_QS, $prev_QE) = ('', '', ''); 
	}
	
	if ( 
	  ( $r1_link_type eq 'p5' and $r1_str eq $r0_str ) 
	  or
	  ( $r1_link_type eq 'p3' and $r1_str ne $r0_str ) 
	  or
	  ( $r1_link_type eq 'incl'                      ) 
	) {
		if ( $r0_str eq 'r' ) {
			&left_closer($r0_QS, $r0_QE, $r1_QS, $r1_QE, $prev_QS, $prev_QE) == 1 and $back_link = $r1_link; 
		} elsif ( $r0_str eq 'f' ) {
			&right_closer($r0_QS, $r0_QE, $r1_QS, $r1_QE, $prev_QS, $prev_QE) == 1 and $back_link = $r1_link; 
		} else {
			die "5: $!\n|$r0_str|\n"; 
		}
	} else {
		return $back_link; 
	}
	
	return $back_link; 
}# End sub closer_r2q_p3()


sub left_closer {
	my ($p0S, $p0E, $p1S, $p1E, $pxS, $pxE) = @_; 
	$p1E <= $p0S or return -1; # This means location ($p1S, $p1E) is not in the left of ($p0S, $p0E) 
	if (defined $pxS and defined $pxE and $pxS ne '') {
		if ( $p1E > $pxE ) {
			return 1; # This means location ($p1S, $p1E) is more closer to ($p0S, $p0E)
		} else {
			return 2; # This means location ($pxS, $pxE) is more closer to ($p0S, $p0E)
		}
	} else {
		return 1; # This means we should use location ($p1S, $p1E); 
	}
}#End sub left_closer() 

sub right_closer {
	my ($p0S, $p0E, $p1S, $p1E, $pxS, $pxE) = @_; 
	$p1S >= $p0E or return -1; # This means location ($p1S, $p1E) is not in the right of ($p0S, $p0E)
	if ( defined $pxS and defined $pxE and $pxS ne '' ) {
		if ( $p1S < $pxS ) {
			return 1; # This means location ($p1S, $#p1E) is closer to ($p0S, $p0E)
		} else {
			return 2; # This means location ($pxS, $pxE) is closer to ($p0S, $p0E)
		}
	} else {
		return 1; # This means we should use location ($p1S, $p1E); 
	}
}#End sub right_closer() 


sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt]", @_); 
}

