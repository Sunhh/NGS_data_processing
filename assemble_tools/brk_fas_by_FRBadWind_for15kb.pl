#!/usr/bin/env perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fastaSunhh; 
my $fasta_obj = fastaSunhh->new(); 
use ConfigSunhh; 
use mathSunhh; 
my $math_obj = mathSunhh->new(); 
use fileSunhh; 
use Getopt::Long; 
my %opts; 

GetOptions(\%opts, 
	'in_fasta:s',  # In_1:  The fasta file to be broken. 
	'in_gapLis:s', # In_2:  $opts{'in_fasta'}.Nlis , the file from deal_fasta.pl -listSite '[nN]+' 
	'in_fwdBad:s', # In_3:  forward_bad_region.wind1kb : key, XXX, WindS, WindE, WindL, BpCnt, FwdBp, RevBp, FwdCnt, RevCnt
	'in_revBad:s', # In_4:  reverse_bad_region.wind1kb : key, XXX, WindS, WindE, WindL, BpCnt, FwdBp, RevBp, FwdCnt, RevCnt
	'in_INSLen:i', # In_5:  Average of insertion size for mapping library. ( 14879 ) 
	'in_maxBadGap:i', # In_6:  Maximum gap size between neighbor bad regions for combining them. (1000)
	'in_maxBad2Gap:i', # In_7:  Maximum distance between directed bad end to a gap for breaking point. (2000)
	'in_kickFlank:i', # In_8:  Flanking size to knock out to uncertain storage. (1000)
	'in_kickCtgRatio:f', # In_9:  Minimum ratio required to set a contig as a false assembly and to replace it with Ns. 
	'in_kickCtgSize:i', # In_10:  Maximum contig length allowed to set a contig as a false assembly and to replace it with Ns. 
	'in_maxGap2BrkPt:i', # In_11:  Maximum length from break point to find a gap and to break that gap. Not used yet. 
	'in_regionToChk:s', # In_12:  A file tell which regions to be checked. 
	'in_brkList:s', # In_13: A file tell how to break the sequence. This will overwrite the -in_fwdBad and -in_revBad options. 
	
	'out_pref:s', # Out_1:  Prefix of output files. Generating {out_pref}.brk.fa and {out_pref}.uncertain.fa
); 

################################################################################
# Setting basic values
################################################################################
$opts{'in_INSLen'} //= 14879 ; # In_5 
$opts{'in_maxBadGap'} //= 1000 ; # In_6
$opts{'in_maxBad2Gap'} //= 2000 ; # In_7
$opts{'in_kickFlank'} //= 1000 ; # In_8
$opts{'in_kickCtgRatio'} //= 0.8 ; # In_9
$opts{'in_kickCtgSize'} //= 10000 ; # In_10
$opts{'in_maxGap2BrkPt'} //= 2000 ; # In_11

$opts{'out_pref'} //= 'out'; # Out_1 

my $is_chk_tse = 0; 
my $ts = 11211000; 
my $te = 11211000; 



################################################################################
# Main body
################################################################################

defined $opts{'in_fasta'} or &usage(); 
defined $opts{'help'} and &usage(); 

my %oriFas_hash = %{ $fasta_obj->save_seq_to_hash( 'faFile'=>$opts{'in_fasta'}, 'has_head'=>1 ) }; 
for my $k1 (keys %oriFas_hash) {
	$oriFas_hash{$k1}{'seq'} =~ s!\s!!g; 
	$oriFas_hash{$k1}{'length'} = length( $oriFas_hash{$k1}{'seq'} ); 
}
my %oriGap_hash; 
if (defined $opts{'in_gapLis'}) {
	my $fh; 
	open $fh,'<',"$opts{'in_gapLis'}" or die; 
	%oriGap_hash = %{ &read_in_gapLis($fh) }; 
	close ($fh); 
} else {
	my $fh; 
	open $fh,'-|',"deal_fasta.pl -listSite '[nN]+' $opts{'in_fasta'}" or die; 
	%oriGap_hash = %{ &read_in_gapLis($fh) }; 
	close ($fh); 
}

my %fwdBad_hash; 
my %revBad_hash; 
my %brk_treatments; 
if (defined $opts{'in_brkList'}) {
	# SeqID   SeqLen  Treat   Num_1   Num_2
	# scaffold5_cov116        1539093 break   2929    3686
	# scaffold5_cov116        1539093 break_flank     171000  171001
	%brk_treatments = %{ &read_in_brkLis( $opts{'in_brkList'} ) }; 
} else {
	%fwdBad_hash = %{ &read_in_badLis( $opts{'in_fwdBad'}, $opts{'in_maxBadGap'} ) }; 
	%revBad_hash = %{ &read_in_badLis( $opts{'in_revBad'}, $opts{'in_maxBadGap'} ) }; 
}

my %out_fh; 
$out_fh{'out_brk'} = &openFH( "$opts{'out_pref'}.brk.fa", '>' ); 
$out_fh{'out_unc'} = &openFH( "$opts{'out_pref'}.uncertain.fa", '>' ); 

print STDOUT join("\t", qw/SeqID SeqLen Treat Num_1 Num_2/)."\n"; 
my $proc_num = 0; 
for my $seqID (sort { $oriFas_hash{$a}{'Order'} <=> $oriFas_hash{$b}{'Order'} } keys %oriFas_hash) {
	my %local; 
	$local{'start'} = 1; 
	$local{'seq_len'} = $oriFas_hash{$seqID}{'length'}; 
	$proc_num ++; 
	&tsmsg("[Msg] Processing $proc_num [$seqID]\n"); 
	my @treatments; 
	if ( defined $opts{'in_brkList'} ) {
		$brk_treatments{ $seqID } //= []; 
		@treatments = @{ $brk_treatments{ $seqID } }; 
	} else {
		@treatments = @{ &find_BrkPt(
		  'seqID' => $seqID, 
		  'oriGap' => \%oriGap_hash, 
		  'fwdBad' => \%fwdBad_hash, 
		  'revBad' => \%revBad_hash, 
		  'local' => \%local, 
		) }; 
	}
	for my $tr (@treatments) {
		print STDOUT join("\t", $seqID, $local{'seq_len'}, @$tr)."\n"; 
	}
	my ($brkSeq_href) = &break_seq( \@treatments, $oriFas_hash{$seqID}{'seq'}, {'in_kickFlank'=>$opts{'in_kickFlank'}, 'in_maxGap2BrkPt'=>$opts{'in_maxGap2BrkPt'}} , $seqID); 
	
	for my $k1 (qw/kept replaced/) {
		defined $brkSeq_href->{$k1} or next; 
		for my $t1 (@{$brkSeq_href->{$k1}}) {
			# ["type:Start-End", seq_string] 
			$t1->[1] =~ s/(.{50})/$1\n/g; chomp($t1->[1]); 
			print {$out_fh{'out_brk'}} ">$seqID:$t1->[0]\n$t1->[1]\n"; 
		}
	}
	defined $brkSeq_href->{'kicked'} or next; 
	for my $t1 (@{$brkSeq_href->{'kicked'}}) {
		# ["type:Start-End", seq_string] 
		$t1->[1] =~ s/(.{50})/$1\n/g; chomp($t1->[1]); 
		print {$out_fh{'out_unc'}} ">$seqID:$t1->[0]\n$t1->[1]\n"; 
	}
}

close ($out_fh{'out_brk'}); 
close ($out_fh{'out_unc'}); 

################################################################################
# Sub-routines 
################################################################################

#	'in_kickFlank:i',    # In_8:   Flanking size to knock out to uncertain storage. (1000)
#	'in_maxGap2BrkPt:i', # In_11:  Maximum length from break point to find a gap and to break that gap. Not used yet. 
# Type: replaceN break_flank break 
#   ['replaceN', $ctg_s, $ctg_e ] 
# Return types : kept replaced kicked 
#        values of each key: [["type:Start-End", seq_string], ...]
#                               type = qw/intact full replaceN break break_flank/ 
sub break_seq {
	my ($brkPt_aref, $oriSeq, $p, $seqID) = @_; 
	my $kickFlank = $p->{'in_kickFlank'}; 
	my $maxGap2BrkPt = $p->{'in_maxGap2BrkPt'}; 
	my $seqLen = length( $oriSeq ); 
	
	my %back; 
	if (scalar( @$brkPt_aref ) == 0) {
		push( @{$back{'kept'} }, ["intact:1-$seqLen", $oriSeq] ); 
		return \%back; 
	}
	# First, do 'replaceN': 
	my @restBrkPt_1; 
	for my $t1 (@$brkPt_aref) {
		my ($type, $s, $e) = @$t1; 
		$type ne 'replaceN' and do { push(@restBrkPt_1, $t1); next; }; 
		push(@{$back{'replaced'}}, ["$type:$s-$e", substr($oriSeq, $s-1, $e-$s+1)]); 
		substr($oriSeq, $s-1, $e-$s+1) = 'N' x ($e-$s+1); 
	}
	if (scalar( @restBrkPt_1 ) == 0) {
		push( @{$back{'kept'}}, ["full:1-$seqLen", $oriSeq] ); 
		return \%back; 
	}
	# Seconcd, break the sequence by order: 
	my $curr_start = 1; 
	my $prev_flank = 0; 
	for my $t1 (@restBrkPt_1) {
		my ($type, $s, $e) = @$t1; 
		$s >= $curr_start or do { $curr_start = $e; next; }; 
		my $prevSeq = substr($oriSeq, $curr_start-1, $s-$curr_start+1); 
		if ( $type eq 'break' ) {
			if ($prev_flank == 1) {
				my $leftFlank_href = &chop_flank_left($prevSeq, $curr_start, $kickFlank); 
				if ($leftFlank_href->{'kicked_seq'} ne '') {
					push(@{$back{'kicked'}}, ["$type:$leftFlank_href->{'kicked_start'}-$leftFlank_href->{'kicked_end'}", $leftFlank_href->{'kicked_seq'}]); 
				}
				if ($leftFlank_href->{'kept_seq'} ne '') {
					$prevSeq = $leftFlank_href->{'kept_seq'}; 
					$curr_start = $leftFlank_href->{'kept_start'}; 
				} else {
					$curr_start = $e; 
					$prev_flank = 0; 
					next; 
				}
			}
			my ($goodPrevSeq, $goodStart, $goodEnd) = &chop_end_Ns($prevSeq, $curr_start); 
			$goodPrevSeq ne '' and push(@{$back{'kept'}}, ["$type:$goodStart-$goodEnd", $goodPrevSeq]); 
			$curr_start = $e; 
			$prev_flank = 0; 
		} elsif ( $type eq 'break_flank' ) {
			if ($prev_flank == 1) {
				my $leftFlank_href = &chop_flank_left($prevSeq, $curr_start, $kickFlank); 
				if ($leftFlank_href->{'kicked_seq'} ne '') {
					push(@{$back{'kicked'}}, ["$type:$leftFlank_href->{'kicked_start'}-$leftFlank_href->{'kicked_end'}", $leftFlank_href->{'kicked_seq'}]); 
				}
				if ($leftFlank_href->{'kept_seq'} ne '') {
					$prevSeq = $leftFlank_href->{'kept_seq'}; 
					$curr_start = $leftFlank_href->{'kept_start'}; 
				} else {
					$curr_start = $e; 
					$prev_flank = 1; 
					next; 
				}
			}
			my $prevSeq_href = &chop_flank_right($prevSeq, $curr_start, $kickFlank); 
			if ($prevSeq_href->{'kept_seq'} ne '') {
				push(@{$back{'kept'}}, ["$type:$prevSeq_href->{'kept_start'}-$prevSeq_href->{'kept_end'}", $prevSeq_href->{'kept_seq'}]); 
			}
			if ($prevSeq_href->{'kicked_seq'} ne '') {
				push(@{$back{'kicked'}}, ["$type:$prevSeq_href->{'kicked_start'}-$prevSeq_href->{'kicked_end'}", $prevSeq_href->{'kicked_seq'}]); 
			}
			$curr_start = $e; 
			$prev_flank = 1; 
		} else {
			&tsmsg("[Err] Bad break_type [$type]\n"); 
		}
	}
	if ( $curr_start <= $seqLen ) {
		my $type = 'AddEnd'; 
		my $prevSeq = substr($oriSeq, $curr_start-1, $seqLen-$curr_start+1); 
		if ($prev_flank == 1) {
			my $leftFlank_href = &chop_flank_left($prevSeq, $curr_start, $kickFlank); 
			if ($leftFlank_href->{'kicked_seq'} ne '') {
				push(@{$back{'kicked'}}, ["$type:$leftFlank_href->{'kicked_start'}-$leftFlank_href->{'kicked_end'}", $leftFlank_href->{'kicked_seq'}]); 
			}
			if ($leftFlank_href->{'kept_seq'} ne '') {
				$prevSeq = $leftFlank_href->{'kept_seq'}; 
				$curr_start = $leftFlank_href->{'kept_start'}; 
			} else {
				$curr_start = $seqLen + 1; 
				$prev_flank = 0; 
			}
		}
		if ( $curr_start <= $seqLen ) {
			my ($goodPrevSeq, $goodStart, $goodEnd) = &chop_end_Ns($prevSeq, $curr_start); 
			$goodPrevSeq ne '' and push(@{$back{'kept'}}, ["$type:$goodStart-$goodEnd", $goodPrevSeq]); 
			$curr_start = $seqLen+1; 
			$prev_flank = 0; 
		}
	}
	
	return \%back; 
}# break_seq 

sub chop_flank_left {
	my ($seq, $start, $kickFlank) = @_; 
	my %back; 
	$back{'start'} = $start; 
	$back{'len'} = length($seq); 
	$back{'end'} = $back{'start'}+$back{'len'}-1; 
	if ( $back{'len'} > $kickFlank ) {
		$back{'kicked_seq'} = substr($seq, 0, $kickFlank); 
		($back{'kicked_seq'}, $back{'kicked_start'}, $back{'kicked_end'}) = &chop_end_Ns( $back{'kicked_seq'}, $back{'start'} ); 
		$back{'kept_seq'}   = substr($seq, $kickFlank, $back{'len'}-$kickFlank); 
		($back{'kept_seq'}, $back{'kept_start'}, $back{'kept_end'}) = &chop_end_Ns( $back{'kept_seq'}, $back{'start'}+$kickFlank-1 ); 
	} else {
		$back{'kept_seq'} = ''; 
		$back{'kept_start'} = -1; 
		$back{'kept_end'} = -1; 
		$back{'kicked_seq'} = $seq; 
		($back{'kicked_seq'}, $back{'kicked_start'}, $back{'kicked_end'}) = &chop_end_Ns( $back{'kicked_seq'}, $back{'start'} ); 
	}
	return \%back; 
}# chop_flank_left() 

sub chop_flank_right {
	my ($seq, $start, $kickFlank) = @_; 
	my %back; 
	$back{'start'} = $start; 
	$back{'len'} = length($seq); 
	$back{'end'} = $back{'start'}+$back{'len'}-1; 
	if ( $back{'len'} > $kickFlank ) {
		$back{'kicked_seq'} = substr($seq, $back{'len'}-$kickFlank, $kickFlank); 
		($back{'kicked_seq'}, $back{'kicked_start'}, $back{'kicked_end'}) = &chop_end_Ns( $back{'kicked_seq'}, $back{'end'}-$kickFlank+1 ); 
		$back{'kept_seq'}   = substr($seq, 0, $back{'len'}-$kickFlank); 
		($back{'kept_seq'}, $back{'kept_start'}, $back{'kept_end'}) = &chop_end_Ns( $back{'kept_seq'}, $back{'start'} ); 
	} else {
		$back{'kept_seq'} = ''; 
		$back{'kept_start'} = -1; 
		$back{'kept_end'} = -1; 
		$back{'kicked_seq'} = $seq; 
		($back{'kicked_seq'}, $back{'kicked_start'}, $back{'kicked_end'}) = &chop_end_Ns( $back{'kicked_seq'}, $back{'end'}-$kickFlank+1 ); 
	}
	return \%back; 
}# chop_flank_right() 

sub chop_end_Ns {
	my ($seq, $start) = @_; 
	my %back; 
	$back{'start'} = $start; 
	$back{'len'} = length($seq); 
	$back{'end'} = $back{'start'}+$back{'len'}-1; 
	if ( $seq =~ s!^([^ATGCatgc]+)!! ) {
		$back{'start'} += length($1); 
	}
	if ( $seq =~ s!([^ATGCatgc]+)$!! ) {
		$back{'end'} -= length($1); 
	}
	if ($back{'start'} > $back{'end'}) {
		return ('', -1, -1); 
	}
	return ($seq, $back{'start'}, $back{'end'}); 
}

sub find_BrkPt {
	my %parm = @_; 
	for my $k1 (qw/seqID oriGap fwdBad revBad local/) {
		defined $parm{$k1} or die "'$k1' is required in &find_BrkPt()\n"; 
	}
	for my $k1 (qw/in_INSLen in_maxBad2Gap in_kickFlank in_kickCtgRatio in_kickCtgSize in_maxGap2BrkPt/) {
		$parm{'local'}{$k1} //= $opts{$k1} ; 
	}
	$parm{'local'}{'start'} //= 1; 
	defined $parm{'local'}{'seq_len'} or die "local{'seq_len'} is required.\n"; 
	
	my @back; 
	my $seqID = $parm{'seqID'}; 
	defined $parm{'fwdBad'}{'loci_array'}{$seqID} or return \@back; 
	FWD_SE: 
	for my $fwd_blk_se (@{$parm{'fwdBad'}{'loci_array'}{$seqID}}) {
		my ($fwd_s, $fwd_e) = @$fwd_blk_se; 
		
		REDO_CURR_FWD: 
		## Enlarge local-start to skip gap. 
		&fix_start(\%parm); 
		## Check step : Re-position fwd_s, fwd_e
		my $record_start = $parm{'local'}{'start'}; 
		$fwd_s < $parm{'local'}{'start'} and $fwd_s = $parm{'local'}{'start'}; 
		$fwd_e < $parm{'local'}{'start'} and next FWD_SE; 
		
		my ($mate_s, $mate_e) = ( $fwd_s, $fwd_e+$parm{'local'}{'in_INSLen'}-1 ); 
		$mate_s > $parm{'local'}{'seq_len'} and next FWD_SE; # This should not happen. 
		$mate_e > $parm{'local'}{'seq_len'} and $mate_e = $parm{'local'}{'seq_len'} ; 
		$mate_s < $parm{'local'}{'start'}+$parm{'local'}{'in_INSLen'}-1 and $mate_s = $parm{'local'}{'start'}+$parm{'local'}{'in_INSLen'}-1; 
		$mate_e < $parm{'local'}{'start'}+$parm{'local'}{'in_INSLen'}-1 and next FWD_SE; # Should not happen. 
		
		if ( my $chk_rule_1 = &pass_rule_1([$mate_s, $mate_e], \%parm) ) {
			# Check if [Fa, Fb+INSLen] overlaps with some [Ra, Rb]s; 
			for my $rev_blk_se (@$chk_rule_1) {
				my ($rev_s, $rev_e) = @$rev_blk_se; 
				$rev_s < $parm{'local'}{'start'}+$parm{'local'}{'in_INSLen'}-1 and $rev_s = $parm{'local'}{'start'}+$parm{'local'}{'in_INSLen'}-1; 
				# Here $rev_e cannot be smaller than $parm{'local'}{'start'}+$parm{'local'}{'in_INSLen'}-1, because $rev_e cannot be smaller than $mate_s; 
				
				if ( my $chk_rule_1_1 = &pass_rule_1_1( [$fwd_s, $fwd_e], [$rev_s, $rev_e], \%parm ) ) {
					if ( $fwd_s <= $rev_s and $fwd_e >= $rev_e ) {
						my $raw_start = $parm{'local'}{'start'}; 
						$parm{'local'}{'start'} = $record_start; 
						if ( my $chk_rule_2_1 = &pass_rule_2_1( [$fwd_s, $fwd_e], \%parm ) ) {
							for my $t1 (@$chk_rule_2_1) {
								my ( $ctg_s, $ctg_e, $ctg_len, $ovl_len ) = @$t1; 
								push(@back, ['replaceN', $ctg_s, $ctg_e]); 
$is_chk_tse and $fwd_s <= $te and $fwd_e >= $ts and &tsmsg("[$fwd_s, $fwd_e] R[$rev_s, $rev_e] rule_1_1_gap\n"); 
								$parm{'local'}{'start'} = &mathSunhh::max( $fwd_e+1, $ctg_e+1 ); 
								goto REDO_CURR_FWD; 
							}
						} else {
							$parm{'local'}{'start'} = $raw_start; 
						}
					}
					push(@back, @$chk_rule_1_1); 
$is_chk_tse and $fwd_s <= $te and $fwd_e >= $ts and &tsmsg("[$fwd_s, $fwd_e] R[$rev_s, $rev_e] rule_1_1\n"); 
					goto REDO_CURR_FWD; 
				} elsif ( my $chk_rule_1_2 = &pass_rule_1_2([$fwd_s, $fwd_e], [$rev_s, $rev_e], \%parm) ) {
					push(@back, @$chk_rule_1_2); 
$is_chk_tse and $fwd_s <= $te and $fwd_e >= $ts and &tsmsg("[$fwd_s, $fwd_e] R[$rev_s, $rev_e] rule_1_2\n"); 
					goto REDO_CURR_FWD; 
				} elsif ( my $chk_rule_1_3 = &pass_rule_1_3([$fwd_s, $fwd_e], [$rev_s, $rev_e], \%parm) ) {
					push(@back, @$chk_rule_1_3); 
$is_chk_tse and $fwd_s <= $te and $fwd_e >= $ts and &tsmsg("[$fwd_s, $fwd_e] R[$rev_s, $rev_e] rule_1_3\n"); 
					goto REDO_CURR_FWD; 
				} else {
					&tsmsg("[Wrn] Unexpected region: [seqID=$seqID, f ($fwd_s, $fwd_e) r ($rev_s, $rev_e)]\n"); 
					next FWD_SE; 
					# die "Why I arrive here!"; 
				}
			}
		} else {
			if ( my $chk_rule_2_1 = &pass_rule_2_1( [$fwd_s, $fwd_e], \%parm ) ) {
				for my $t1 (@$chk_rule_2_1) {
					my ( $ctg_s, $ctg_e, $ctg_len, $ovl_len ) = @$t1; 
					# [ $ctg_s, $ctg_e, $ctg_len, $ovl_len ] 
					push(@back, ['replaceN', $ctg_s, $ctg_e ]); # Sometimes $ctg_s <= 'start', which means in the new sequence, this is a break instead of replaceN. 
					$parm{'local'}{'start'} = &mathSunhh::max( $fwd_e+1, $ctg_e+1 ); 
				}
$is_chk_tse and $fwd_s <= $te and $fwd_e >= $ts and &tsmsg("[$fwd_s, $fwd_e] rule_2_1\n"); 
				goto REDO_CURR_FWD; 
			} else {
				# $parm{'local'}{'start'} = $fwd_e+1; # Because the OriSeq is not changed, so the usable start position should remain unchaged. 
				next FWD_SE; 
			}
		}# 
	}# for:FWD_SE 
	
	
	return \@back; 
}# find_BrkPt() 

sub fix_start {
	my ($p) = @_; 
	my $seqID = $p->{'seqID'}; 
	defined $p->{'oriGap'}{'loci_array'}{$seqID} or return undef(); 
	for my $gSE (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
		$gSE->[0] > $p->{'local'}{'start'} and last; 
		$gSE->[1] < $p->{'local'}{'start'} and next; 
		$gSE->[0] <= $p->{'local'}{'start'} and $gSE->[1] >= $p->{'local'}{'start'} and $p->{'local'}{'start'} = $gSE->[1]+1; 
	}
	return undef(); 
}# fix_start() 

sub pass_rule_1_3 {
	my ($se1_aref, $se2_aref, $p) = @_; 
	my ($fa, $fb) = @$se1_aref; 
	my ($ra, $rb) = @$se2_aref; 
	# $se1_aref = [Fa, Fb] ; $se2_aref = [Ra, Rb] ; 
	( $fa <= $ra and $ra <= $fb ) or return undef(); 
	my $seqID = $p->{'seqID'}; 
	my @back; 
	
	if ($ra == $fb) {
		push(@back, ['break_flank', $ra-1, $fa]); 
		$p->{'local'}{'start'} = $fa + $p->{'local'}{'in_kickFlank'}; 
	} else {
		push(@back, ['break_flank', $ra-1, $ra]); 
		push(@back, ['break_flank', $fb, $fb+1]); 
		$p->{'local'}{'start'} = $fb+1 + $p->{'local'}{'in_kickFlank'}; 
	}
	
	return \@back; 
}# pass_rule_1_3() 

sub pass_rule_2_1 {
	my ($se1_aref, $p) = @_; 
	my $seqID = $p->{'seqID'}; 
	my $seqLen = $p->{'local'}{'seq_len'}; 
	my $maxCtgLen = $p->{'local'}{'in_kickCtgSize'}; 
	my $minCtgRatio = $p->{'local'}{'in_kickCtgRatio'}; 
	my @ctg_se; 
	my $se2_aaref = $p->{'oriGap'}{'loci_array'}{$seqID}; 
	if (defined $se2_aaref and @$se2_aaref > 0) {
		for (my $i=0; $i<$#{$se2_aaref}; $i++) {
			my $ctg_s = $se2_aaref->[$i][1]+1; 
			my $ctg_e = $se2_aaref->[$i+1][0]-1; 
			$ctg_e < $se1_aref->[0] and next; 
			$ctg_s > $se1_aref->[1] and last; 
			$ctg_s < $p->{'local'}{'start'} and $ctg_s = $p->{'local'}{'start'}; 
			my $ctg_len = $ctg_e-$ctg_s+1; 
			my ($ovl_len) = $math_obj->ovl_region( $se1_aref->[0], $se1_aref->[1], $ctg_s, $ctg_e ); 
			push(@ctg_se, [ $ctg_s, $ctg_e, $ctg_len, $ovl_len ]); 
		}
		my ($ctg_s, $ctg_e) = ($se2_aaref->[-1][1]+1, $seqLen); 
		$ctg_s < $p->{'local'}{'start'} and $ctg_s = $p->{'local'}{'start'}; 
		if ( $ctg_e >= $se1_aref->[0] and $ctg_s <= $se1_aref->[1] ) {
			my $ctg_len = $ctg_e-$ctg_s+1; 
			my ($ovl_len) = $math_obj->ovl_region( $se1_aref->[0], $se1_aref->[1], $ctg_s, $ctg_e ); 
			push(@ctg_se, [ $ctg_s, $ctg_e, $ctg_len, $ovl_len ]); 
		}
	} else {
		@ctg_se = ( [$p->{'local'}{'start'}, $seqLen, $seqLen-$p->{'local'}{'start'}+1, $se1_aref->[1]-$se1_aref->[0]+1] ); 
	}
	my @back; 
	for my $t1 (@ctg_se) {
		$t1->[2] > $maxCtgLen and next; 
		($t1->[2]-1000) * $minCtgRatio > $t1->[3] and next; 
		push(@back, [@$t1]); 
	}
	scalar(@back) > 0 or return undef(); 
	return \@back; 
}# pass_rule_2_1() 

sub pass_rule_1 {
	my ($se1_aref, $p) = @_; 
	my $seqID = $p->{'seqID'}; 
	
	my @back; 
	for my $t1 (@{$p->{'revBad'}{'loci_array'}{$seqID}}) {
		$t1->[0] > $se1_aref->[1] and last; 
		$t1->[1] < $se1_aref->[0] and next; 
		# ( defined $t1->[0] and defined $t1->[1] ) or die "t1=@$t1\n"; 
		push(@back, [@$t1]); 
	}
	scalar(@back) > 0 or return undef(); 
	if (scalar(@back) > 1 and $back[0][1]-$back[0][0]+1 <= 1001) {
		for my $t1 (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
			$t1->[0] > $back[0][1] and last; 
			$t1->[1] < $back[0][0] and next; 
			&mathSunhh::min($t1->[1], $back[0][1]) - &mathSunhh::max($t1->[0], $back[0][0]) + 1 >= 500 and do { shift(@back); last; }; 
		}
	}
	return \@back; 
}# pass_rule_1() 

# 				if ( my $chk_rule_1_1 = &pass_rule_1_1( [$fwd_s, $fwd_e], [$rev_s, $rev_e], \%parm ) ) {
sub pass_rule_1_1 {
	my ($se1_aref, $se2_aref, $p) = @_; 
	my ($fa, $fb) = @$se1_aref; 
	my ($ra, $rb) = @$se2_aref; 
	# $se1_aref = [Fa, Fb] ; $se2_aref = [Ra, Rb] ; 
	my $maxBad2Gap = $p->{'local'}{'in_maxBad2Gap'}; 
	my $seqID = $p->{'seqID'}; 
	
	my @back; 
	defined $p->{'oriGap'}{'loci_array'}{$seqID} or return undef(); 
	
	# The first choice is to break from the one close to rev_bad_blk_start; 
	# The second choice is to break from the one close to fwd_bad_blk_end; 
	if ( $fb < $rb ) {
		my $wind_len = 1000; 
		for my $t1 (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
			my ($gs, $ge) = @{$t1}; 
			$ge < $p->{'local'}{'start'} and next; 
			$ge < $fa and $ge < $ra-$maxBad2Gap and next; 
			$gs > $fb+$maxBad2Gap and $gs > $rb+$maxBad2Gap and last; 
			$gs < $p->{'local'}{'start'} and die "Why arrives here! seqID=$seqID gap[$gs, $ge], Fab[$fa, $fb], Rab[$ra, $rb], start=[$p->{'local'}{'start'}]\n"; 
			
			if ( $ra-$fa+1 <= $p->{'local'}{'in_INSLen'}+$wind_len and $ge >= $ra-$maxBad2Gap and $ge < $ra+$maxBad2Gap) {
				push(@back, ['break', $gs-1, $ge+1]); 
				$p->{'local'}{'start'} = $ge+1; 
				return \@back; 
				last; 
			}
		}
		for my $t1 (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
			my ($gs, $ge) = @{$t1}; 
			$ge < $p->{'local'}{'start'} and next; 
			$ge < $fa and $ge < $ra-$maxBad2Gap and next; 
			$gs > $fb+$maxBad2Gap and $gs > $rb+$maxBad2Gap and last; 
			$gs < $p->{'local'}{'start'} and die "Why arrives here! seqID=$seqID gap[$gs, $ge], Fab[$fa, $fb], Rab[$ra, $rb], start=[$p->{'local'}{'start'}]\n"; 
			
			if ( $fb-$fa+1 <= $p->{'local'}{'in_INSLen'}+$wind_len and $gs <= $fb+$maxBad2Gap and $gs > $fb-$maxBad2Gap) { 
				push(@back, ['break', $gs-1, $ge+1]); 
				$p->{'local'}{'start'} = $ge+1; 
				return \@back; 
				last; 
			}
		}
	}
	
	
	for my $t1 (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
		my ($gs, $ge) = @{$t1}; 
		$ge < $p->{'local'}{'start'} and next; 
		$ge < $fa and $ge < $ra-$maxBad2Gap and next; 
		$gs > $fb+$maxBad2Gap and $gs > $rb+$maxBad2Gap and last; 
		$gs < $p->{'local'}{'start'} and die "Why arrives here! seqID=$seqID gap[$gs, $ge], Fab[$fa, $fb], Rab[$ra, $rb], start=[$p->{'local'}{'start'}]\n"; 
		
		if ( $fb < $ra and $gs > $ra+$maxBad2Gap ) {
			# Rule_1_1_1: 
			my $chk_inner_rule_1_2 = &pass_rule_1_2( [$fa, $fb], [$ra, $rb], $p ); 
			push(@back, @$chk_inner_rule_1_2); 
			last; 
		} elsif ( $fa <= $ra and $ra <= $fb ) {
			# Rule_1_1_2: 
			# Break at (Ra-1, Ra) and (Fb, Fb+1)
			if ( $ra == $fb ) {
				push(@back, ['break_flank', $ra-1, $fb+1]); # It impossible that $ra==start, because I recompute Ra according to local-start; 
			} else {
				push(@back, ['break_flank', $ra-1, $ra]); # It impossible that $ra==start, because I recompute Ra according to local-start; 
				push(@back, ['break_flank', $fb, $fb+1]); # 
			}
			$p->{'local'}{'start'} = $fb+1+$p->{'local'}{'in_kickFlank'}; 
			last; 
		} else {
			# Rule_1_1_3: 
			push(@back, ['break', $gs-1, $ge+1]); 
			$p->{'local'}{'start'} = $ge+1; 
			last; 
		}
	}
	
	scalar(@back) > 0 or return undef(); 
	return \@back; 
}# pass_rule_1_1() 
#	for my $k1 (qw/seqID oriGap fwdBad revBad local/) {
#		defined $parm{$k1} or die "'$k1' is required in &find_BrkPt()\n"; 
#	}

# pass_rule_1_2( [$fa, $fb], [$ra, $rb], $p ); 
sub pass_rule_1_2 {
	my ($se1_aref, $se2_aref, $p) = @_; 
	my ($fa, $fb) = @$se1_aref; 
	my ($ra, $rb) = @$se2_aref; 
	$fb < $ra or return undef(); 
	# $se1_aref = [Fa, Fb] ; $se2_aref = [Ra, Rb] ; 
	my $seqID = $p->{'seqID'}; 
	my @back; 
	
	# First I need to find $FnextS and $RprevE; 
	# For $FnextS
	my $FnextS = &get_FnextS([$fa, $fb], $p); 
	my $RprevE = &get_RprevE([$ra, $rb], $p); 
	
	if ( $FnextS-$fb >= $ra-$RprevE ) {
		# [Fb, GoodE] is longer. 
		push(@back, ['break_flank', $fb, $fb+1]); 
		$p->{'local'}{'start'} = $fb+1+$p->{'local'}{'in_kickFlank'}; 
	} else {
		# [GoodS, Fa] is longer. 
		push(@back, ['break_flank', $ra-1, $ra]); 
		$p->{'local'}{'start'} = $ra+$p->{'local'}{'in_kickFlank'}; 
	}
	
	return \@back; 
}# pass_rule_1_2() 
#	for my $k1 (qw/seqID oriGap fwdBad revBad local/) {
#		defined $parm{$k1} or die "'$k1' is required in &find_BrkPt()\n"; 
#	}
#	for my $k1 (qw/in_INSLen in_maxBad2Gap in_kickFlank in_kickCtgRatio in_kickCtgSize in_maxGap2BrkPt/) {
#		$parm{'local'}{$k1} //= $opts{$k1} ; 
#	}

# I don't consider overlapping here!!! 
sub get_FnextS {
	my ($fAB_aref, $p) = @_; 
	my ($fa, $fb) = @$fAB_aref; 
	my $seqID = $p->{'seqID'}; 
	my $fnextS; 
	for my $fSE (@{$p->{'fwdBad'}{'loci_array'}{$seqID}}) {
		my ($fs, $fe) = @$fSE; 
		$fs <= $fb and next; 
		$fnextS = $fs-1; 
		last; 
	}
	$fnextS //= $p->{'local'}{'seq_len'}; 
	for my $gSE (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
		my ($gs, $ge) = @$gSE; 
		$gs <= $fb and next; 
		$fnextS > $gs-1 and $fnextS = $gs-1; 
		last; 
	}
	
	return $fnextS; 
}# get_FnextS() 

# I don't consider overlapping here!!! 
sub get_RprevE {
	my ($rAB_aref, $p) = @_; 
	my ($ra, $rb) = @$rAB_aref; 
	my $seqID = $p->{'seqID'}; 
	my $rprevE; 
	for my $rSE (@{$p->{'revBad'}{'loci_array'}{$seqID}}) {
		my ($rs, $re) = @$rSE; 
		$re >= $ra and last; 
		$rprevE = $re+1; 
	}
	$rprevE //= $p->{'local'}{'start'}; 
	$rprevE < $p->{'local'}{'start'} and $rprevE = $p->{'local'}{'start'}; 
	for my $gSE (@{$p->{'oriGap'}{'loci_array'}{$seqID}}) {
		my ($gs, $ge) = @$gSE; 
		$ge >= $ra and last; 
		$rprevE < $ge+1 and $rprevE = $ge+1; 
	}
	
	return $rprevE; 
}# get_RprevE() 

sub read_in_brkLis {
	my ($fn) = @_; 
	my %back; 
	open F,'<',"$fn" or die "Failed to open $fn. $!\n"; 
	# SeqID   SeqLen  Treat   Num_1   Num_2
	# scaffold5_cov116        1539093 break   2929    3686
	# scaffold5_cov116        1539093 break_flank     171000  171001
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$. == 1 and $ta[0] =~ m/^SeqID$/i and next; 
		my ($id, $seqLen, $treat, $n1, $n2) = @ta; 
		push(@{$back{$id}}, [$treat, $n1, $n2]); 
	}
	close F; 
	return \%back; 
}# End read_in_brkLis() 

sub read_in_badLis {
	my ($fn, $dist2jn) = @_; 
	$dist2jn //= $opts{'in_maxBadGap'}; 
	my %back; 
	open F,'<', "$fn" or die; 
	# [zf25@cbsulm02 03.compare_errors]$ head tag_fwd.out
	# ChromID BpCnt   WindS   WindE   WindL   BpCnt.1 FwdBp   RevBp   FwdCnt  RevCnt
	# scaffold116_cov122      1000    41001   42000   1000    1000    483     0       0       19
	# scaffold116_cov122      1000    42001   43000   1000    1000    1000    853     0       24
	# scaffold116_cov122      1000    43001   44000   1000    1000    1000    1000    0       24
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$. == 1 and ( $ta[0] =~ m/^(ChromID|ChrID)$/i ) and next; 
		my ($id, $s, $e) = @ta[0,2,3]; 
		push(@{$back{'loci_array'}{$id}}, [$s, $e]); 
	}
	close F; 
	for my $k1 (keys %{$back{'loci_array'}}) {
		$back{'loci_array'}{$k1} = $math_obj->mergeLocBlk( $back{'loci_array'}{$k1}, 'dist2join'=>$dist2jn+1 ); 
		$back{'loci_to_idx'}{$k1} = &_sort_index_loci_array( $back{'loci_array'}{$k1} ); 
	}
	
	return \%back; 
}# read_in_badLis() 

sub read_in_gapLis {
	my ($fh) = @_; 
	my %back; 
	# Key     Length  MatchStart      MatchEnd        MatchLen
	# scaffold2_cov61 1823    1363    1630    268
	# scaffold5_cov116        1539093 2930    3685    756
	# scaffold5_cov116        1539093 5846    7254    1409
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$. == 1 and $ta[0] eq 'Key' and next; 
		my ($id, $len, $s, $e) = @ta[0,1,2,3]; 
		$s > $e and die "s=$s , e=$e , s should not be larger than e.\n"; 
		push(@{$back{'loci_array'}{$id}}, [$s, $e]); 
		# $back{'loci_to_idx'}{$id}{"$s-$e"} = $#{ $back{'loci_array'}{$id} }; 
		$back{'seq_len'}{$id} //= $len; 
	}
	for my $id (keys %{$back{'loci_array'}}) {
		$back{'loci_to_idx'}{$id} = &_sort_index_loci_array( $back{'loci_array'}{$id} ); 
	}
	return \%back; 
}# read_in_gapLis() 

# Return : {"$s-$e"} => $index_in_array 
#          Also sort and change input loci_array. 
sub _sort_index_loci_array {
	my ($aref_0) = @_; 
	@$aref_0 = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$aref_0; 
	my %back; 
	for (my $i=0; $i<@$aref_0; $i++) {
		my $tk = "$aref_0->[$i][0]-$aref_0->[$i][1]"; 
		defined $back{$tk} and die "Repeated tk=[$tk]\n";
		$back{$tk} = $i; 
	}
	return \%back; 
}

sub usage {
	print <<HH; 
################################################################################
perl $0 ... 
	'in_fasta:s',  # In_1:  The fasta file to be broken. 
	'in_gapLis:s', # In_2:  \$opts{'in_fasta'}.Nlis , the file from deal_fasta.pl -listSite '[nN]+' 
	'in_fwdBad:s', # In_3:  forward_bad_region.wind1kb : key, XXX, WindS, WindE, WindL, BpCnt, FwdBp, RevBp, FwdCnt, RevCnt
	'in_revBad:s', # In_4:  reverse_bad_region.wind1kb : key, XXX, WindS, WindE, WindL, BpCnt, FwdBp, RevBp, FwdCnt, RevCnt
	'in_INSLen:i', # In_5:  Average of insertion size for mapping library. ( 14879 ) 
	'in_maxBadGap:i', # In_6:  Maximum gap size between neighbor bad regions for combining them. (1000)
	'in_maxBad2Gap:i', # In_7:  Maximum distance between directed bad end to a gap for breaking point. (2000)
	'in_kickFlank:i', # In_8:  Flanking size to knock out to uncertain storage. (1000)
	'in_kickCtgRatio:f', # In_9:  Minimum ratio required to set a contig as a false assembly and to replace it with Ns. 
	'in_kickCtgSize:i', # In_10:  Maximum contig length allowed to set a contig as a false assembly and to replace it with Ns. 
	'in_maxGap2BrkPt:i', # In_11:  Maximum length from break point to find a gap and to break that gap. 
	'in_brkList:s', # In_13: A file tell how to break the sequence. This will overwrite the -in_fwdBad and -in_revBad options.
	
	'out_pref:s', # Out_1:  Prefix of output files. Generating {out_pref}.brk.fa and {out_pref}.uncertain.fa
################################################################################
HH
	exit 1; 
}# usage() 
