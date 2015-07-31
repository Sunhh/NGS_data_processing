#!/usr/bin/perl
# 20150729 I still need to check for repeats' alignments. 
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"minIdentity:f", # Default 0, range [0-100]
	"minIdentityWhenOvl:f", # Default 0, range [0-100], but I might use 88.97 or something. 
	"repeat_dist:f", 
	"useRawID!", 
); 

-t and !@ARGV and die "perl $0 filtered_sorted.Spec97_to_ScfPI.bn6\n"; 

$opts{'minIdentity'}        //= 0; 
$opts{'minIdentityWhenOvl'} //= 0; 
$opts{'repeat_dist'}        //= -1; 

my %p1_to_p2; 
my %p1_to_p2_repeat; 
my %len1; 
my %len2; 

while (<>) {
	$. % 100e3 == 1 and &tsmsg("[Msg] Processing $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	my $identity = $ta[2]; 
	$identity >= $opts{'minIdentity'} or next; 
	my ( $id1, $id2, $l1, $l2 ) = @ta[0,1,12,13]; 
	my ( $s1,$e1, $s2, $e2 ) = @ta[6,7,8,9]; 
	my $str = $ta[14]; 
	$len1{$id1} //= $l1; 
	$len2{$id2} //= $l2; 
#	my ($tag1, $scfID_1, $scfS_1, $scfE_1) = &parse_specID( $id1 ); 
#	&add_pos_toHash(
#	  \%p1_to_p2, 
#	  [$scfID_1, $scfS_1, $scfE_1], 
#	  [$id1, $s1, $e1, $id2, $s2, $e2], 
#	  [$identity, $opts{'minIdentityWhenOvl'}], 
#	  \%p1_to_p2_repeat
#	); 
	&add_pos_toHash(
	  \%p1_to_p2, 
	  [$id1, $s1, $e1], 
	  [$id1, $s1, $e1, $id2, $s2, $e2], 
	  [$identity, $opts{'minIdentityWhenOvl'}], 
	  \%p1_to_p2_repeat, 
	  $opts{'repeat_dist'}
	); 
}

&tsmsg("[Rec] Output table.\n"); 
#                                    0            1          2       3          4        5   6        ?
print STDOUT join("\t", qw/FromScfID FromScfStart FromScfEnd ToScfID ToScfStart ToScfEnd Str Identity Tag/)."\n"; 
for my $scfID1 (sort keys %p1_to_p2) {
	for my $ar1 (@{$p1_to_p2{$scfID1}}) {
		my ( $s1, $e1, $id2, $s2, $e2, $str, $ident ) = @$ar1; 
		#ar1-[0], [1], [2] , [3], [4], [5],  [6]
		if ($opts{'useRawID'}) {
			unless ( defined $p1_to_p2_repeat{$scfID1} ) {
				my ( $origID, $origS, $origE ) = &SpecLoc_to_OrigLoc($scfID1, $s1, $e1); 
				print STDOUT join("\t", $origID, $origS, $origE, $id2, $s2, $e2, $str, $ident, 'Unique')."\n"; 
				next; 
			}
			$p1_to_p2_repeat{$scfID1} //= []; 
			my ($restSE, $ovlSE, $idx_start) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( $p1_to_p2_repeat{$scfID1}, [$s1, $e1] ); 
			for my $ar2 ( @$restSE ) {
				my ($restS, $restE) = @$ar2; 
				#ar2-[0],   [1]
				my ($new_S2, $new_E2) = &mathSunhh::map_loc_byReference(
				  [ $s1, $e1, $s2, $e2 ], 
				  [ $restS, $restE ], 
				  $str
				); 
				my ( $origID, $origS, $origE ) = &SpecLoc_to_OrigLoc($scfID1, $restS, $restE); 
				print STDOUT join("\t", $origID, $origS, $origE, $id2, $new_S2, $new_E2, $str, $ident, "Unique")."\n"; 
			}
			for my $ar3 ( @$ovlSE ) {
				my ($ovlS, $ovlE) = @$ar3; 
				#ar3-[0],  [1]
				my ( $origID, $origS, $origE ) = &SpecLoc_to_OrigLoc($scfID1, $ovlS, $ovlE); 
				print STDOUT join("\t", $origID, $origS, $origE, $id2, '.', '.', $str, $ident, "Repeat")."\n"; 
			}
		} else {
			unless ( defined $p1_to_p2_repeat{$scfID1} ) {
				print STDOUT join("\t", $scfID1, $s1, $e1, $id2, $s2, $e2, $str, $ident, 'Unique')."\n"; 
				next; 
			}
			$p1_to_p2_repeat{$scfID1} //= []; 
			my ($restSE, $ovlSE, $idx_start) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( $p1_to_p2_repeat{$scfID1}, [$s1, $e1], [6] ); 
			for my $ar2 ( @$restSE ) {
				my ($restS, $restE) = @$ar2; 
				#ar2-[0],   [1]
				my ($new_S2, $new_E2) = &mathSunhh::map_loc_byReference(
				  [ $s1, $e1, $s2, $e2 ], 
				  [ $restS, $restE ], 
				  $str
				); 
				print STDOUT join("\t", $scfID1, $restS, $restE, $id2, $new_S2, $new_E2, $str, $ident, "Unique")."\n"; 
			}
			for my $ar3 ( @$ovlSE ) {
				my ($ovlS, $ovlE) = @$ar3; 
				#ar3-[0],  [1]
				print STDOUT join("\t", $scfID1, $ovlS, $ovlE, $id2, '.', '.', $str, $ident, "Repeat")."\n"; 
			}
		}
	}
}

&tsmsg("[Rec] Finished. $0\n"); 

# add_pos_toHash(\%p1_to_p2, [$scfID_1, $scfS_1, $scfE_1], [$id1, $s1, $e1, $id2, $s2, $e2], [$identity, $opts{'minIdentityWhenOvl'}]); 
#	&add_pos_toHash(
#	  \%p1_to_p2, 
#	  [$id1, $s1, $e1], 
#	  [$id1, $s1, $e1, $id2, $s2, $e2], 
#	  [$identity, $opts{'minIdentityWhenOvl'}], 
#	  \%p1_to_p2_repeat, 
#	  $opts{'repeat_dist'}
#	); 
sub add_pos_toHash {
	my ($href, $aref_1, $aref_2, $aref_3, $href_repeat, $dist_ident) = @_; 
	$aref_3 //= [-1, -1]; 
	$aref_3->[0] //= -1; # Current identity%
	$aref_3->[1] //= -1;   # Cutoff for second identity%
	$href_repeat //= {};  # The hash storing repeat (multiple aligned) regions. 
	$dist_ident  //= -1;  # Cutoff for recognizing repeats. Distance of current and previous identity%. 
	my ( $scfID_1, $scfS_1, $scfE_1 ) = @$aref_1; 
	my ( $id1, $s1, $e1, $id2, $s2, $e2 ) = @$aref_2; 
	my $str = ( $s2 <= $e2 ) ? '+' : '-' ; 
	$s2 > $e2 and ($s2, $e2) = ($e2, $s2); 
	unless ( defined $href->{$scfID_1} ) {
		push(@{$href->{$scfID_1}}, [$scfS_1, $scfE_1, $id2, $s2, $e2, $str, $aref_3->[0]]); 
		return; 
	}
	#my $new_S1 = $s1 + $scfS_1 - 1; 
	#my $new_E1 = $e1 + $scfS_1 - 1; 
	my ($restSE, $ovlSE, $idx_start) = &mathSunhh::sep_loc2_by_loc1_singleLoc2( $href->{$scfID_1}, [$scfS_1, $scfE_1], [6] ); 
	## Firstly, check for repeat region. 
	if ( $dist_ident >= 0 and $aref_3->[0] > 0 and @$ovlSE > 0 ) {
		for my $ar2 (@$ovlSE) {
			$aref_3->[0] >= $ar2->[2]-$dist_ident or next; 
			&add_pos_toHash(
			  $href_repeat, 
			  [$scfID_1, $ar2->[0], $ar2->[1]], 
			  [$scfID_1, $ar2->[0], $ar2->[1], $id2, $s2, $e2], 
			  [200, 0], 
			  {}, 
			  -1
			); 
		}
	}
	## Then check for specific region. 
	@$restSE == 0 and return; 
	if ( $aref_3->[1] > 0 or @$restSE > 1 or $restSE->[0] != $scfS_1 or $restSE->[1] != $scfE_1 ) {
		$aref_3->[0] >= $aref_3->[1] or return; 
	}
	for my $ar1 (@$restSE) {
		my ($new_S2, $new_E2) = &mathSunhh::map_loc_byReference(
		  [$scfS_1,$scfE_1,$s2,$e2], 
		  [$ar1->[0], $ar1->[1]], 
		  $str
		); 
		# push(@{$href->{$scfID_1}}, [$ar1->[0], $ar1->[1], $id2, $new_S2, $new_E2, $str]); 
		&mathSunhh::insert_SEloc_toSEArray( $href->{$scfID_1}, [ $ar1->[0], $ar1->[1], $id2, $new_S2, $new_E2, $str, $aref_3->[0] ], $idx_start ); 
	}
	@{$href->{$scfID_1}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$href->{$scfID_1}}; 
	return ; 
}# add_pos_toHash() 

sub parse_specID {
	my $in = shift; 
	$in =~ m/^(Ref[^\s_]+)_(\S+)_(\d+)_(\d+)$/ or &stopErr("[Err] Unknown SpecID [$in]\n"); 
	my ($refTag, $scfID, $scfS, $scfE) = ($1, $2, $3, $4); 
	return($refTag, $scfID, $scfS, $scfE); 
}
# my ( $origID, $origS, $origE ) = &SpecLoc_to_OrigLoc($scfID1, $ovlS, $ovlE); 
# Return   : ( $origID, $origS, $origE )
sub SpecLoc_to_OrigLoc {
	my ($specID, $specS, $specE) = @_; 
	my ( $origID, $origS, $origE ); 
	my ( $refTag, $refID, $refS, $refE ) = &parse_specID( $specID ); 
	$specS //= 1; 
	$specE //= ($refE-$refS+1); 
	$origID = $refID; 
	$origS = $refS + $specS - 1; 
	$origE = $refS + $specE - 1; 
	return ( $origID, $origS, $origE ); 
}# SpecLoc_to_OrigLoc() 

