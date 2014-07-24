#!/usr/bin/perl
use strict;
use warnings; 
use Getopt::Long;

my %opts;
GetOptions(\%opts,
	"maxBlkGap:i", 
	"maxScfGap:i", 
	"inCtgFa:s", 
	"specs:s", 
	"noReplace!", 
	"outScfFas:s", 
	"outLnkTbl:s", 
	"outScfLnk:s", 
	"scfPref:s", 
	"help!", 
); 

# 0                  1               2       3       4               5               6             7               8               9               10              11              12      13              14
# Ref1ID             Ref2ID          Ref1Len Ref2Len Ref1_blk_normS  Ref1_blk_normE  Ref1_blk_fr   Ref2_blk_normS  Ref2_blk_normE  Ref2_blk_fr     Gap_bwt_Blocks  QryID           QryLen  QryBlk1         QryBlk2
# Ref.C6112509__0.0  Ref.scaffold48  500     119264  1               500             r             4260            46522           f               11190           Qry.scaffold    61662   r:7560-8059     f:19250-60747
# Ref.C6112547__0.0  Ref.scaffold37  500     81960   1               459             r             8607            10795           f               1895            Qry.scaffold73  58746   f:56954-57412   r:52820-55058
# Ref.C6112593__0.0  Ref.scaffold12  500     1168    1               451             r             1               1168            r               5584            Qry.scaffold14  186114  f:6945-7395     f:307-1360

-t and !@ARGV and &usage(); 
$opts{help} and &usage(); 

#####################################################
# Help documents. 
sub usage {
	print STDOUT <<HELP; 
######################################################################################################
perl $0 in.link

-help    
-maxBlkGap    [15e3] Maximum gap length accepted between two nearby blocks. 
-maxScfGap    [15e3] Maximum gap length accepted between two nearby scaffolds estimated. 
-outLnkTbl    [\\*STDOUT] Write accepted links to be used for scaffolding. 
-inCtgFa      [AllNamedCtg.fas] Should have "spec." tag consistant with -specs
-specs        [Ref,Qry] Prefix in the beginning of contig names from two sets. 
-noReplace    [False] When given, will not try to find better sequence to replace the spec1 set. 

-outScfLnk    [\\*STDOUT] Write scaffold links 
-outScfFas    [File name] Write scaffold sequences to this file if given. 
-scfPref      [NULL] Prefix of output scaffolds name. 

######################################################################################################
HELP
	exit 1; 
}


#####################################################
# Basic VARs 
# 
my $glob_noReplace = 0; 
$opts{noReplace} and $glob_noReplace = 1; 
defined $opts{scfPref} or $opts{scfPref} = ''; 
defined $opts{maxBlkGap} or $opts{maxBlkGap} = 15000; 
defined $opts{maxScfGap} or $opts{maxScfGap} = 15000; 
defined $opts{specs} or $opts{specs} = "Ref,Qry"; 
my ($spec_1, $spec_2) = split(/,/, $opts{specs}); 
$spec_1 =~ s/^\s+//; $spec_1 =~ s/\s+$//; 
$spec_2 =~ s/^\s+//; $spec_2 =~ s/\s+$//; 
my $has_fa = 0; 
my (%infor_1, %infor_2); 
my %seqLength; 
if (defined $opts{inCtgFa}) {
	&tsmsg("[Rec] Reading -inCtgFa $opts{inCtgFa}\n"); 
	open FA,'<',"$opts{inCtgFa}" or die; 
	$has_fa = 1; 
	my $tk; 
	my $tlist; 
	while (<FA>) {
		if (m/^\s*\>/) {
			m/^\s*\>(\S+)/ or die "FA_line:$_\n"; 
			$tk = $1; 
			if ( $tk =~ m/^$spec_1\.\S+/ ) {
				$tlist = \%{$infor_1{seq}}; 
			} elsif ( $tk =~ m/^$spec_2\.\S+/ ) {
				$tlist = \%{$infor_2{seq}}; 
			} else {
				die "FA_name:$_\n"; 
			}
		} else {
			$tlist->{$tk} .= $_; 
		}
	}
	close FA; 
	&tsmsg("[Rec] Preparing information for contigs.\n"); 
	for ( keys %{$infor_1{seq}} ) {
		$infor_1{seq}{$_}     =~ s/\s//g; 
		$infor_1{seq}{$_}     = uc( $infor_1{seq}{$_} ); 
		$infor_1{atgcLis}{$_} = &listATGC( \$infor_1{seq}{$_} ); 
		$infor_1{len}{$_}     = length( $infor_1{seq}{$_} ); 
		$seqLength{$_}        = $infor_1{len}{$_}; 
	}
	for ( keys %{$infor_2{seq}} ) {
		$infor_2{seq}{$_}     =~ s/\s//g; 
		$infor_2{seq}{$_}     = uc( $infor_2{seq}{$_} ); 
		$infor_2{atgcLis}{$_} = &listATGC( \$infor_2{seq}{$_} ); 
		$infor_2{len}{$_}     = length( $infor_2{seq}{$_} ); 
		# $seqLength{$_}        = $infor_2{len}{$_}; 
	}
}

my $oScfFasFh = undef(); 
if (defined $opts{outScfFas}) {
	my $tfh; 
	open $tfh, '>', "$opts{outScfFas}" or die; 
	$oScfFasFh = $tfh; 
}
my $oLnkTblFh = \*STDOUT; 
if (defined $opts{outLnkTbl}) {
	my $tfh; open $tfh, '>', "$opts{outLnkTbl}" or die; 
	$oLnkTblFh = $tfh; 
}
my $oScfLnkFh = \*STDOUT; 
if (defined $opts{outScfLnk}) {
	my $tfh; 
	open $tfh, '>', "$opts{outScfLnk}" or die; 
	$oScfLnkFh = $tfh; 
}


my %blk_link; 
#   %{$r_link} : {scfID_1}{fr_1}{scfID_2} = [ fr_2, [scfLen_1,scfS_1,scfE_1], [scfLen_2,scfS_2,scfE_2], [Blk_GapLen, Scf_GapLen], [QryID,QryLen, QryFR_1,QryS_1,QryE_1, QryFR_2,QryS_2,QryE_2] ]
my %mess_link; # Paired links supporting both "A-> B->" and "A-> <-B". 




#####################################################
# Read in in.link file 
#
&tsmsg("[Rec] Reading link file.\n"); 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[0] eq 'Ref1ID' and next; 
	&get_link(\%blk_link, \%mess_link, @ta); 
	defined $seqLength{$ta[0]} and $seqLength{$ta[0]} != $ta[2] and die "Diff_len: $seqLength{$ta[0]} for $ta[0] : $_\n"; 
	defined $seqLength{$ta[1]} and $seqLength{$ta[1]} != $ta[3] and die "Diff_len: $seqLength{$ta[1]} for $ta[1] : $_\n"; 
	defined $seqLength{$ta[0]} or $seqLength{$ta[0]} = $ta[2];
	defined $seqLength{$ta[1]} or $seqLength{$ta[1]} = $ta[3];  
}
&tsmsg("[Rec] link file has been read in.\n"); 

#####################################################
# Check for conflict paired-links, and remove them. 
# %mess_link || %mult_link 
#
{
my %mult_link = %{ &multiLink(\%blk_link) }; # Paired links supporting at least "A-> B" and "A-> C"
&tsmsg("[Rec] Drop [" . &delete_link2(\%blk_link, \%mess_link) . "] mess links and [" . &delete_link2(\%blk_link, \%mult_link) . "] multi-links.\n"); 
my %mmmm_link = %{ &multiLink(\%blk_link) }; # Paired links supporting at least "A-> B" and "A-> C"
&tsmsg("[Rec] rest " . scalar(keys %mmmm_link) . "\n"); 
}

#####################################################
# Output links 
# &printLink(\%blk_link); 

#####################################################
# Build new scaffolds
# 
&tsmsg("[Rec] Building scaffolds.\n"); 
my ($scfLnk, $scfSeq) = &build_scaffold( \%blk_link, \%seqLength, \%infor_1, \%infor_2, $glob_noReplace); 
## Output scf links 
if ( defined $oScfLnkFh and scalar(@$scfLnk) > 0 ) {
	&tsmsg("[Rec] Writing scaffold link record.\n"); 
	my $nn = int( log( scalar(@$scfLnk) )/log(10) ) + 1; 
	for (@$scfLnk) {
		my $sc_name = sprintf("%s%0${nn}d", $opts{scfPref}, $_->[0]); 
		print {$oScfLnkFh} ">$sc_name $_->[1]\n$_->[2]\n"; 
	}
}
## Output scf sequences. 
if ( defined $oScfFasFh and $has_fa and scalar(@$scfSeq) > 0 ) {
	&tsmsg("[Rec] Writing scaffold sequences.\n"); 
	my $nn = int( log( scalar(@$scfSeq) )/log(10) ) + 1; 
	for ( @$scfSeq ) {
		my $sc_name = sprintf("%s%0${nn}d", $opts{scfPref}, $_->[0]); 
		my $sc_len = length($_->[1]); 
		print {$oScfFasFh} ">$sc_name $sc_len\n$_->[1]\n"; 
	}
}


#####################################################
# All done! 
&tsmsg("[Rec] All done!\n"); 


#####################################################
# Sub-functions 
#

# Input  : ( $infor_1{atgcLis}{scfID_Ref}, $infor_1{len}{scfID_Ref}, start, end )
# Return : ratio_of_ATGC_in_region 
sub ratioATGC {
	my ($ar, $seqLen, $start, $end) = @_; 
	defined $start or $start = 1; 
	defined $end   or $end = $seqLen; 
	$start > $end and return -1; 

	my $sum = 0; 
	for (@$ar) {
		$_->[0] > $end and last; 
		$_->[1] < $start and next; 
		my ($ts, $te) = @$_; 
		$ts < $start and $ts = $start; 
		$te > $end and $te = $end; 
		$sum += ($te-$ts+1); 
	}

	my $ratio = sprintf("%0.4f", $sum/($end-$start+1)); 
	return $ratio; 
}#End sub ratioATGC() 

# Return : [ [matchS1, matchE1], [matchS2, matchE2], ... ]
sub listATGC ($) {
	my $sr = shift; 
	my @back; 
	
	pos($$sr) = 0; 
	while ( $$sr =~ m/\G(?:.*?)([ATGCatgc]+)/gs ) {
		push(@back, [ $-[1]+1 , $+[1] ]); 
		pos($$sr) = $+[1]; 
	}
	pos($$sr) = 0; 

	return \@back; 
}#End sub listATGC() 

sub build_scaffold {
	my ($pair, $ctgLen, $infor1, $infor2, $no_replace) = @_; 
	my ($sc_ct); 
	my $seen_start = {}; # A hash reference. 
	my @back_ScfLnk; 
	my @back_ScfSeq; 
	
	# Build scaffolds from seed contig
	SEED: 
	for my $scfID_1 (sort { $ctgLen->{$b} <=> $ctgLen->{$a} } keys %$ctgLen) {

		defined $seen_start->{$scfID_1} and next SEED; 
		
		# From $scfID_1 to 3p direction; 
		my @chain_right = ( [ $scfID_1, 'f', $ctgLen->{$scfID_1} ] ); 
		@chain_right = &computeLayout( "R", \@chain_right, [ $scfID_1, 'f', $ctgLen->{$scfID_1} ], $pair, $ctgLen, $seen_start ); 

		# From $scfID_1 to 5p direction; 
		my @chain_left; 
		@chain_left = &computeLayout( "L", \@chain_left, [ $scfID_1, 'r', $ctgLen->{$scfID_1} ], $pair, $ctgLen, $seen_start ); 

		my @chain_total = (@chain_left, @chain_right); 

		# Output scaffold 
		$sc_ct ++; 
		$sc_ct % 1000 == 1 and &tsmsg("[Msg] Have $sc_ct scaffolds.\n"); 
		# scalar(@chain_total) > 1 or next SEED; 
		push(@back_ScfLnk, &printScaffLink("${sc_ct}", \@chain_total) ); 
		if ( defined $oScfFasFh and $has_fa ) {
			my $sc_seq = &getScaffSeq( \@chain_total, $infor1, $infor2, $no_replace ); 
			push(@back_ScfSeq, [ $sc_ct, $sc_seq ]); 
		}

	}#End (SEED:) for my $scfID_1 

	return (\@back_ScfLnk, \@back_ScfSeq); 
}#End sub build_scaffold() 


# Input  : ( 
#            \%infor_1_Ref, 
#            \%infor_2_Qry, 
#            $prev_Scf_E, 
#            [ #0_prev_Scf_ID, #1_prev_Scf_FR, #2_prev_Scf_length ], 
#            [ #0_scfGap_1to2, #1_blkGap_1to2, 
#              #2_prev_Blk_S,  #3_prev_Blk_E, #4_foll_Blk_S, #5_foll_Blk_E, 
#              #6_QryID,       #7QryLen, 
#              #8_QryFR_prev,  #9_QryS_prev,  #10_QryE_prev, 
#              #11_QryFR_foll, #12_QryS_foll, #13_QryE_foll 
#            ], 
#            $no_replace
#          )
# Return : ( $seq_to_add, $End_position_of_seq )
sub seqWithBlk {
	my ($infor1, $infor2, $prev_Scf_E, $ta1, $ta2, $no_replace) = @_; 
	defined $no_replace or $no_replace = 0; 
	my @ta1 = @$ta1; 
	my @ta2 = @$ta2; 
	my $back_seq = ''; 

	if ( $ta1[1] eq 'f' ) {
		if ( $ta2[3] < $prev_Scf_E ) {
			# Bad case 1a : Impossible. 
			die "Bad case 1: $!\n"; 
		} elsif ( $ta2[3] == $prev_Scf_E ) {
			return ('', $prev_Scf_E); 
		} else {
			if ( $ta2[2] <= $prev_Scf_E ) {
				$back_seq = &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], $prev_Scf_E+1, $ta2[3] ); 
				return ( $back_seq, $ta2[3] ); 
			} else {
				$back_seq  = &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], $prev_Scf_E+1, $ta2[2]-1 ); 
				my ($seq2add, $tE) = &betterSeq( 
				  [ $ta1[0], @ta2[2,3], $ta1[1] ], 
				  [ @ta2[6,  9,10,      8]      ], 
				  $infor1, $infor2, $no_replace
				); 
				$back_seq .= $seq2add; 
				$tE == $ta2[3] or die "$!\n"; 
				return ( $back_seq, $tE ); 
			}
		}
	} elsif ( $ta1[1] eq 'r' ) {
		if ( $ta2[2] > $prev_Scf_E ) {
			# Bad case 1b : Impossible. 
			die "Bad case 1b: $!\n $ta2[2] > $prev_Scf_E\n@ta1\n@ta2\n"; 
		} elsif ( $ta2[2] == $prev_Scf_E ) {
			return ('', $prev_Scf_E); 
		} else {
			if ( $ta2[3] >= $prev_Scf_E ) {
				$back_seq = &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], $ta2[2], $prev_Scf_E-1 ); 
				return( $back_seq, $ta2[2] ); 
			} else {
				$back_seq  = &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], $ta2[3]+1, $prev_Scf_E-1 ); 
				my $tmp_str = $ta2[8]; $tmp_str =~ tr/fr/rf/; 
				my ($seq2add, $tE) = &betterSeq(
				  [ $ta1[0], @ta2[2,3], $ta1[1]   ],
				  [ $ta2[6], @ta2[9,10], $tmp_str ], 
				  $infor1, $infor2, $no_replace
				); 
				$back_seq .= $seq2add; 
				$tE == $ta2[2] or die "$!\n"; 
				return( $back_seq, $tE ); 
			}
		}
	} else {
		die "$!\n"; 
	}
	die "$!\n"; 
}#End sub seqWithBlk

# Input  : ( \@scf_chain_total, \%infor_1_spec1_Ref, \%infor_2_spec2_Qry, "0/1"_If_not_replace_by_better_seq )
# Return : $scafSequence 
sub getScaffSeq {
	my ($chain, $infor1, $infor2, $no_replace) = @_; 
	my $newseq = ''; 
	defined $no_replace or $no_replace = 0; 
	my $prev_Scf_E = -1; 
	for (my $i=0; $i<@$chain; $i+=2) {

		my @ta1 = @{$chain->[$i]};   # ( #0_prev_Scf_ID, #1_prev_Scf_FR, #2_prev_Scf_length ) 

		### For the first scaffold, the prev_Scf_E needs to be changed. 
		if ( $prev_Scf_E == -1 ) {
			$i == 0 or die "Bad case 2:\n"; 
			$ta1[1] eq 'f' and $prev_Scf_E = 0; 
			$ta1[1] eq 'r' and $prev_Scf_E = $infor1->{len}{$ta1[0]}+1; 
		}

		### When we meet the last scaffold, there is no following block left. 
		if ( $i == $#{$chain} ) {
			if ( $ta1[1] eq 'f' ) {
				$newseq .= &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], $prev_Scf_E+1, $infor1->{len}{$ta1[0]} ); 
			} elsif ( $ta1[1] eq 'r' ) {
				$newseq .= &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], 1,             $prev_Scf_E-1           ); 
			} else {
				die "$!\n"; 
			}
			last; 
		}
		$i == $#{$chain} and do { $newseq .= &subSeq( \$infor1->{seq}{$ta1[0]}, $ta1[1], $prev_Scf_E+1 ); die "Bad here.\n"; last; }; 

		my @ta2 = @{$chain->[$i+1]}; # ( #0_scfGap_1to2, #1_blkGap_1to2, 
		                           #   #2_prev_Blk_S,  #3_prev_Blk_E, #4_foll_Blk_S, #5_foll_Blk_E, 
		                           #   #6_QryID,       #7QryLen, 
		                           #   #8_QryFR_prev,  #9_QryS_prev,  #10_QryE_prev, 
		                           #   #11_QryFR_foll, #12_QryS_foll, #13_QryE_foll )

		my @ta3 = @{$chain->[$i+2]}; # ( #0_foll_Scf_ID, #1_foll_Scf_FR, #2_foll_Scf_length ) 

		# Add sequence from prev_Scf_E to the end of the previous block. 
		my ($add_seq_with_blk, $prev_blkE) = &seqWithBlk( $infor1, $infor2, $prev_Scf_E, \@ta1, \@ta2, $no_replace ); 
		$newseq .= $add_seq_with_blk; 

		# Add gap sequence from outer of the previous block to the outer of the following block. 
		my ($add_GapSeq, $gapE) = &betterGapSeq(
		  [ $ta1[0], @ta2[2,3],  $ta1[1], $ta3[0], @ta2[4,5],   $ta3[1], $ta2[0] ], 
		  [ $ta2[6], @ta2[9,10], $ta2[8],          @ta2[12,13], $ta2[11]         ], 
		  $infor1, $infor2, $no_replace
		); 
		$newseq .= $add_GapSeq; 

		# Add block sequence from the following contig. 
		my $tmp_str = $ta2[11]; 
		$ta3[1] eq 'r' and $tmp_str =~ tr/fr/rf/; 
		my ($add_blkSeq, $foll_blkE) = &betterSeq(
		  [ $ta3[0], @ta2[4,5], $ta3[1] ], 
		  [ @ta2[6,  12,13],    $tmp_str], 
		  $infor1, $infor2, $no_replace
		); 
		$newseq .= $add_blkSeq; 
		$prev_Scf_E = $foll_blkE; 


	}#End for (my $i=0; $i<@$chain; $i+=2) 
	
	return $newseq; 
}#End sub getScaffSeq() 

# Input  : ( [id1, s1, e1, str1], [id2, s2, e2, str2], infor1, infor2, "0/1"_If_not_replace ) : betterSeq() 
#          str1/str2 should be same to the direction laying on the scaffold chain. 
#          So str2 should not be directly retrived from chain. 
# Return : ( $seq_to_add, $last_position_of_seq_along_scaf_chain ) 
sub betterSeq {
	my ($reg1, $reg2, $infor1, $infor2, $no_replace) = @_; 
	defined $no_replace or $no_replace = 0; 
	my ($id1, $s1, $e1, $str1) = @$reg1; 
	my ($id2, $s2, $e2, $str2) = @$reg2; 
	( $e1 >= $s1 and $e2 >= $s2 ) or die "asfd\n"; 
	my $back_seq; 
	my $tE; 
	$str1 eq 'f' and $tE = $e1; 
	$str1 eq 'r' and $tE = $s1; 

	if ( $no_replace ) {
		$back_seq = &subSeq( \$infor1->{seq}{$id1}, $str1, $s1, $e1 ); 
		return( $back_seq, $tE ); 
	} 

	# Check which one is better. 
	## Rule 1. The one with larger ATGC% is better. 
	my $atgcR_1 = &ratioATGC( $infor1->{atgcLis}{$id1}, $infor1->{len}{$id1}, $s1, $e1 ); 
	my $atgcR_2 = &ratioATGC( $infor2->{atgcLis}{$id2}, $infor2->{len}{$id2}, $s2, $e2 ); 
	if ( $atgcR_1 >= $atgcR_2 ) {
		$back_seq = &subSeq( \$infor1->{seq}{$id1}, $str1, $s1, $e1 ); 
	} elsif ( $atgcR_1 < $atgcR_2 ) {
		$back_seq = &subSeq( \$infor2->{seq}{$id2}, $str2, $s2, $e2); 
	} else {
		;# Not here. 
	}

	return ($back_seq, $tE); 
}# sub betterSeq() 

# Input  : ( 
#           [$id1a, $s1a_blkS, $e1a_blkE, $str1a, $id1b, $s1b, $e1b, $str1b, $scfGap], 
#           [$id2,  $s2a_blkS, $e2a_blkE, $str2a,        $s2b, $e2b, $str2b         ], 
#           \%infor1, \%infor2, 
#           "0/1"_If_not_replace ) : betterGapSeq() 
# Return : ( $seq_to_add, $End_position_of_seq )
sub betterGapSeq {
	my ( $reg1, $reg2, $infor1, $infor2, $no_replace ) = @_; 
	defined $no_replace or $no_replace = 0; 
	my ( $id1a, $s1a, $e1a, $str1a, $id1b, $s1b, $e1b, $str1b, $scfGap ) = @$reg1; 
	my ( $id2,  $s2a, $e2a, $str2a,        $s2b, $e2b, $str2b          ) = @$reg2; 
	( $e1a >= $s1a+1 and $e1b >= $s1b+1 ) or die "afsd\n"; 
	my $back_seq; 
	my $tE; 

	# Compute Gap location of previous and following contigs. 
	my ( $gapS1a, $gapE1a, $gapS1b, $gapE1b, $gapSeq ); 
	my ( $gapS2, $gapE2, $gapStr2a, $gapStr2b ); 
	$gapStr2a = $str2a; $gapStr2b = $str2b; 
	if ( $str1a eq 'f' )    { $gapS1a=$e1a+1; $gapE1a=$infor1->{len}{$id1a};                          } 
	elsif ( $str1a eq 'r' ) { $gapS1a=1;      $gapE1a=$s1a-1;                $gapStr2a =~ tr/fr/rf/;  } 
	else                    { die "$!\n"; }
	$gapS1a <= $gapE1a + 1 or die "$!\n"; 
	if ( $str1b eq 'f' )    { $gapS1b=1;      $gapE1b=$s1b-1;                                        $tE = $gapE1b; } 
	elsif ( $str1b eq 'r' ) { $gapS1b=$e1b+1; $gapE1b=$infor1->{len}{$id1b}; $gapStr2b =~ tr/fr/rf/; $tE = $gapS1b; } 
	else                    { die "$!\n"; }
	$gapS1b <= $gapE1b+1 or die "$!\n"; 

	$gapStr2a eq $gapStr2b or die "gap strand not same!\n@$reg1\n@$reg2\n"; 
	if    ( $gapStr2a eq 'f' ) { $gapS2 = $e2a+1; $gapE2 = $s2b-1; } 
	elsif ( $gapStr2a eq 'r' ) { $gapS2 = $e2b+1; $gapE2 = $s2a-1; }
	else                       { die "$!\n"; }
	$gapS2 <= $gapE2+1 or die "$!\n"; 

	my $gapNSeq = ( $scfGap < 0 ) ? "N" : "N" x $scfGap ; 
	my $gapNLen = ( $scfGap < 0 ) ? 1   : $scfGap ; 
	
	if ($no_replace) {
		$back_seq  = &subSeq( \$infor1->{seq}{$id1a}, $str1a, $gapS1a, $gapE1a ); 
		$back_seq .= $gapNSeq; 
		$back_seq .= &subSeq( \$infor1->{seq}{$id1b}, $str1b, $gapS1b, $gapE1b ); 
		return ($back_seq, $tE); 
	}

	# Check which one is better. 
	my $atgcR_1a = &ratioATGC( $infor1->{atgcLis}{$id1a}, $infor1->{len}{$id1a}, $gapS1a, $gapE1a ); 
	my $atgcR_1b = &ratioATGC( $infor1->{atgcLis}{$id1b}, $infor1->{len}{$id1b}, $gapS1b, $gapE1b ); 
	my $atgcR_2  = &ratioATGC( $infor2->{atgcLis}{$id2},  $infor2->{len}{$id2},  $gapS2,  $gapE2  ); 
	my ($gapLen1a, $gapLen1b, $gapLen2, $gapLen1); 
	$gapLen1a = $gapE1a-$gapS1a+1; 
	$gapLen1b = $gapE1b-$gapS1b+1; 
	$gapLen2  = $gapE2-$gapS2+1; 
	$atgcR_1a == -1 and $gapLen1a = 0; 
	$atgcR_1b == -1 and $gapLen1b = 0; 
	$atgcR_2  == -1 and $gapLen2  = 0; 
	$gapLen1  = $gapLen1a + $gapNLen + $gapLen1b; 
	my $atgcR_1  = ( $gapLen1 == 0 ) ? -1 : ( $atgcR_1a * $gapLen1a + $atgcR_1b * $gapLen1b ) / $gapLen1 ; 
	## Rule 1. No gap will be the best. 
	$gapLen2 == 0 and return ('', $tE); 
	## Rule 2. The one with larger ATGC% is better. 
	if ( $atgcR_2 >= $atgcR_1 ) {
		$back_seq = &subSeq( \$infor2->{seq}{$id2}, $gapStr2a, $gapS2, $gapE2 ); 
	} else {
		$back_seq  = &subSeq( \$infor1->{seq}{$id1a}, $str1a, $gapS1a, $gapE1a ); 
		$back_seq .= $gapNSeq; 
		$back_seq .= &subSeq( \$infor1->{seq}{$id1b}, $str1b, $gapS1b, $gapE1b ); 
	}
	return ($back_seq, $tE); 

}#End sub betterGapSeq() 

# Input  : (\$seq, $str, $start, $end) : subSeq() 
sub subSeq {
	my ($seqR, $str, $s, $e) = @_; 
	my $subs; 

	$s < 1 and return ''; 
	if (defined $e) {
		$e >= $s or return ''; 
		$subs = substr($$seqR, $s-1, $e-$s+1); 
	} else {
		$subs = substr($$seqR, $s-1); 
	}
	if ( $str eq 'f' ) {
		; 
	} elsif ( $str eq 'r' ) {
		$subs = reverse($subs); 
		$subs =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
	} else {
		die "str[$str] not known.\n"; 
	}
	return $subs; 
}#End subSeq()

sub printScaffLink {
	my ($sc_name, $chain) = @_; 
	my $back_txt = ''; 
	my $sum = 0; 
	for my $tr (@$chain) {
		@$tr == 3 or next; 
		$sum += $tr->[2]; 
	}
	# $back_txt .= ">$sc_name $sum\n"; 
	for my $tr (@$chain) {
		$back_txt .= join(' ', @$tr)."\n"; 
	}
	return [$sc_name, $sum, $back_txt]; 
}#End sub printScaffLink 

sub computeLayout {
	my ( $ext, $chain, $init_tig, $pair, $ctgLen, $seen_start ) = @_; 
	my $extension = 1; 
	my $curr_tig = [@$init_tig]; 

	EXTENSION: 
	while ( $extension ) {
		# Exit if $curr_tig has been included in some scaffold. 
		defined $seen_start->{ $curr_tig->[0] } and do { $extension = 0; last EXTENSION; }; 

		# Set $curr_tig->[0] as included if the init_tig is [ scfID_1, 'r' ] and . 
		# Must be careful when searhcing for the next_tig, because when the init_tig is [ scfID_1, 'f' ], $seen_start->{$scfID_1} is always "undef()"; 
		( $curr_tig->[0] eq $init_tig->[0] and $init_tig->[1] eq 'f' ) or $seen_start->{ $curr_tig->[0] } ++; 

		# Exit if there is no following contig (next_tig) defined in %$pair. 
		( defined $pair->{ $curr_tig->[0] } and defined $pair->{ $curr_tig->[0] }{ $curr_tig->[1] } ) or do { $extension = 0; last EXTENSION; }; 

		
		my $list = $pair->{ $curr_tig->[0] }{ $curr_tig->[1] }; 
		my $next_tig = []; 
		$next_tig->[0] = (keys %$list)[0];                  # scfID_2 : Following contig. 
		$next_tig->[1] = $list->{ $next_tig->[0] }; 
		# my $next_array = $list->{ $next_tig->[0] }; 
		# $next_tig->[1] = $next_array->[0];                  # scfFR-2 
		# $next_tig->[2] = $next_array->[2][0];   # scfLen_2
		# $next_tig->[3] = $next_array->[3][1];   # scfGap_1to2 after scfID_1 
		# $next_tig->[4] = $next_array->[3][0];   # blkGap_1to2 after scfID_1 

		# Some checks to make sure the following link is good. 
		#### Should not be the original input contig ($init_tig->[0]) 
		$next_tig->[0] eq $init_tig->[0] and do { $extension = 0; last EXTENSION; }; 
		#### Should not be used in other scaffolds.  
		defined $seen_start->{ $next_tig->[0] } and do { $extension = 0; last EXTENSION; }; 
		#### Gap length ? Or other rules. 

		# push the following contig (next_tig) to the chain. 
		$curr_tig = &getChain($ext, $chain, $next_tig); 

	}#End while ($extension)

	return @$chain; 
}#End sub computeLayout() 

# Return  : [ $next_tig_id, $next_tig_fr, $next_tig_ll/length ]
# Changed : 
#    @$chain : ( 
#               [ prev_Scf_ID, prev_Scf_FR, prev_Scf_length ], 
#               [ scfGap_1to2, blkGap_1to2, prev_Blk_S, prev_Blk_E, foll_Blk_S, foll_Blk_E, QryID,QryLen, QryFR_prev,QryS_prev,QryE_prev, QryFR_foll,QryS_foll,QryE_foll ], 
#               [ foll_Scf_ID, foll_Scf_FR, foll_Scf_length ], 
#               ....
#              )
sub getChain {
	my ($ext, $chain, $next_tig) = @_; 
	my $next_tig_id = $next_tig->[0]; 
	my @tb = @{$next_tig->[1]}; # Value of %{$r_link}
	my $next_tig_fr = $tb[0]; 
	my $next_tig_ll = $tb[2][0]; 

	if ( $ext eq 'R' ) {
		# 3p direction. 
		push(@$chain, 
		  [ $tb[3][1], $tb[3][0], @{$tb[1]}[1,2], @{$tb[2]}[1,2], @{$tb[4]},  ], 
		  [ $next_tig_id, $next_tig_fr, $next_tig_ll ]
		); 
		# push(@$chain, [@ta[3,4,0,1,2]]); 
	} elsif ( $ext eq 'L' ) {
		# 5p direction. 
		my $tmp_fr = $next_tig_fr; 
		$tmp_fr =~ tr/fr/rf/; 
		unshift(@$chain, 
		  [ $next_tig_id, $tmp_fr, $next_tig_ll ], 
		  [ $tb[3][1], $tb[3][0], @{$tb[2]}[1,2], @{$tb[1]}[1,2], @{$tb[4]}[0,1, 5,6,7, 2,3,4] ]
		); 
		# unshift(@$chain, [@ta]); 
	} else {
		die "why here?\n"; 
	}
	return [ $next_tig_id, $next_tig_fr, $next_tig_ll ]; 
}#End sub getChain() 


#   %{$r_link} : {scfID_1}{fr_1}{scfID_2} = [ fr_2, [scfLen_1,scfS_1,scfE_1], [scfLen_2,scfS_2,scfE_2], [Blk_GapLen, Scf_GapLen], [QryID,QryLen, QryFR_1,QryS_1,QryE_1, QryFR_2,QryS_2,QryE_2] ]
sub printLink {
	my $r = shift; 
	print $oLnkTblFh join("\t", qw/BLK_GapLen SCF_GapLen BLK_ID_1 BLKID_2 LINK Qry_ID_1 Qry_ID_2/)."\n"; 
	for my $scfID_1 (sort keys %$r) {
		for my $fr_1 (sort keys %{$r->{$scfID_1}}) {
			for my $scfID_2 (sort keys %{$r->{$scfID_1}{$fr_1}}) {
				my $rA = $r->{$scfID_1}{$fr_1}{$scfID_2}; 
				my $fr_2 = $rA->[0]; 

				my $blkID_1 = "BLK:$scfID_1:$rA->[1][0]:$fr_1:$rA->[1][1]:$rA->[1][2]"; 
				my $blkID_2 = "BLK:$scfID_2:$rA->[2][0]:$fr_2:$rA->[2][1]:$rA->[2][2]"; 
				my $blkGap  = "BLK_Gap:$rA->[3][0]"; 
				my $scfGap  = "SCF_Gap:$rA->[3][1]"; 
				my $qryID_1 = "REF:$rA->[4][0]:$rA->[4][1]:$rA->[4][2]:$rA->[4][3]:$rA->[4][4]"; 
				my $qryID_2 = "REF:$rA->[4][0]:$rA->[4][1]:$rA->[4][5]:$rA->[4][6]:$rA->[4][7]"; 
				print $oLnkTblFh join("\t", $blkGap, $scfGap, $blkID_1, $blkID_2, "LINK", $qryID_1, $qryID_2)."\n"; 
			}
		}
	}
	return 0; 
}#End sub printLink() 

# (\%src_link_to_delete, \%index_to_delete)
# Delete hash_{key1}{key2} from %src_link_to_delete which exist in %index_to_delete . 
sub delete_link2 {
	my ($srcL, $idxL) = @_; 
	my $cnt = 0; 
	for my $k1 ( keys %$idxL ) {
		for my $k2 ( keys %{$idxL->{$k1}} ) {
			exists $srcL->{$k1}{$k2} or next; 
			delete $srcL->{$k1}{$k2}; 
			$cnt ++; 
		}
	}
	return $cnt; 
}#End sub delete_link2() 




sub multiLink {
	my $r_link = shift; 
	my %bad_scf; 
	for my $scf1 (keys %$r_link) {
		for my $fr ( keys %{$r_link->{$scf1}} ) {
			scalar(keys %{$r_link->{$scf1}{$fr}}) > 1 and do { $bad_scf{$scf1}{$fr} = 1; }; 
		}
	}
	return \%bad_scf; 
}


# Output: 
#   %{$r_link} : {scfID_1}{fr_1}{scfID_2} = [ fr_2, [scfLen_1,scfS_1,scfE_1], [scfLen_2,scfS_2,scfE_2], [Blk_GapLen, Scf_GapLen], [QryID,QryLen, QryFR_1,QryS_1,QryE_1, QryFR_2,QryS_2,QryE_2] ]
sub get_link {
	my ($r_link, $mess_link, @ta) = @_; 
	# Record in forward direction. 
	my $tA = &getLinkArray(\@ta); 
	defined $tA or return 0; # This link cannot pass GapLen limit. 

	if ( defined $r_link->{$ta[0]} and defined $r_link->{$ta[0]}{$ta[6]} and defined $r_link->{$ta[0]}{$ta[6]}{$ta[1]}) {
		$tA->[0] eq $ta[9] or do { &tsmsg("[Wrn] Multiple conflict links for {$ta[0]}{$ta[6]}{$ta[1]}.\n"); $mess_link->{$ta[0]}{$ta[6]}{$ta[1]} = 1; }; 
	} else {
		$r_link->{$ta[0]}{$ta[6]}{$ta[1]} = $tA; 
	}

	# Record in reverse direction. 
	my @tb = @ta[1,0, 3,2, 7,8,9, 4,5,6, 10,11,12, 14,13]; 
	$tb[6] =~ tr/fr/rf/; $tb[9] =~ tr/fr/rf/; 

	if ( defined $r_link->{$tb[0]} and defined $r_link->{$tb[0]}{$tb[6]} and defined $r_link->{$tb[0]}{$tb[6]}{$tb[1]} ) {
		$r_link->{$tb[0]}{$tb[6]}{$tb[1]}[0] eq $tb[9] or do { &tsmsg("[Wrn] Multiple conflict links for {$tb[0]}{$tb[6]}{$tb[1]}.\n"); $mess_link->{$tb[0]}{$tb[6]}{$tb[1]} = 1; }; 
	} else {
		$r_link->{$tb[0]}{$tb[6]}{$tb[1]} = &getLinkArray(\@tb); 
	}

	return 0; 
}#End sub get_link() 

#   %{$r_link} : {scfID_1}{fr_1}{scfID_2} = [ fr_2, [scfLen_1,scfS_1,scfE_1], [scfLen_2,scfS_2,scfE_2], [Blk_GapLen, Scf_GapLen], [QryID,QryLen, QryFR_1,QryS_1,QryE_1, QryFR_2,QryS_2,QryE_2] ]
sub getLinkArray {
	my $ar=shift; 
	my @back; 
	$back[0] = $ar->[9]; 
	my @ta1 = ( $ar->[13] =~ m/^([fr]):(\d+)\-(\d+)$/ ) or die "Parse: $ar->[13]\n"; # (QryFR,QryS,QryE)
	my @ta2 = ( $ar->[14] =~ m/^([fr]):(\d+)\-(\d+)$/ ) or die "Parse: $ar->[14]\n"; 
	$back[1] = [ $ar->[2], $ar->[4], $ar->[5]]; 
	$back[2] = [ $ar->[3], $ar->[7], $ar->[8] ]; 
	my $rs1 = ( $ar->[6] eq 'f' ) ? $ar->[2]-$ar->[5] : $ar->[4]-1       ; # Rest length of scaffold 1
	my $rs2 = ( $ar->[9] eq 'f' ) ? $ar->[7]-1        : $ar->[3]-$ar->[8]; # Rest length of scaffold 2
	$back[3] = [ $ar->[10], $ar->[10]-$rs1-$rs2  ]; # Block_Gap , Scaff_Gap
	$back[3][0] >= $opts{maxBlkGap} and return undef(); 
	$back[3][1] >= $opts{maxScfGap} and return undef(); 
	$back[4] = [ $ar->[11], $ar->[12], @ta1[0..2], @ta2[0..2] ]; 
	return \@back; 
}#End sub getLinkArray() 

sub tsmsg {
	my $tt = scalar(localtime()); 
	print STDERR join('', "[$tt]", @_); 
}
