#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use gffSunhh; 
use fastaSunhh; 
use fileSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	# Input files 
	"inGff:s", "scfFa:s", "seqInGff!", 
	
	# Output files 
	"out:s", 
	"silent!", 
	
	# Actions 
	"seqret!", 
	
	"compare2gffC:s", 
	 "sameAll!", 
	 "sameIntron!", "sameSingleExon!", 
	 "rmOvlap!", "rmOvlapLen:i", "rmOvlapRatio1:f", "rmOvlapRatio2:f", "rmOvlapType:s", "rmOvlapStrand:s", 
	
	"ovlapLongest!", "ovlapLength:i", "ovlapRatio:f", "ovlapStrand:s", "ovlapFeatType:s", 
	"islandGene:i", "islandFeatType:s", "islandStrand:s", 
	"gffret:s", "idType:s", 
	"listTopID!", 
	
	"addFaToGff!", 
	
	# Filter options. 
	"extractFeat:s", 
	"help!", 
); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 -inGff in.gff 
# 
# -help             Give help information. 
# 
# Input options: 
# -inGff            [] input.gff3 . requried. 
# -scfFa            [] fasta sequences 
# -seqInGff         [Boolean] Get sequence from gff file if given. 
# 
# Output options: 
# -out              [\*STDOUT] Output file. 
# -sortGffBy        [lineNum] Could be 'linenum/str_posi/position' . 
#                     lineNum  : Output offspring gff lines within each topID in the input order. 
#                     str_posi : Sort offspring lines by 'scfID' and 'strand' and 'position'. 
#                     position : Sort offspring lines by 'scfID' and 'position' only. 
# 
# Action options: 
#----------------------------------------------------------------------------------------------------
# -addFaToGff       [Boolean] Add -scfFa to the tail of -inGff . 
#----------------------------------------------------------------------------------------------------
# -seqret           [Boolean] Retrieve fasta sequence. Not used yet. 
#
#----------------------------------------------------------------------------------------------------
# -gffret           [input_ID_list] I use the first column. 
#   -idType         [mRNA] Could be 'gene/exon/CDS/match/match_part'
#----------------------------------------------------------------------------------------------------
# -listTopID        [Boolean] 
#   -idType         [Same to -gff_top_hier if not given]. 
#----------------------------------------------------------------------------------------------------
# -compare2gffC     [] Name of gff3 file (gffC) to compare with. 
#
#   -sameAll          [Boolean]
#   
#   -sameIntron       [Boolean] Keep multi-exon gene models with same intron composition. 
#                       For single-exon gene models, only require the 1st gene is part of the 2nd (gffC). 
#     -sameSingleExon   [Boolean] Require single-exon gene model exactly the same, too. 
#                         If this is not given, the inGff will be output if it is included by gffC. 
#  
#   -rmOvlap          [Boolean] Remove 
#     -rmOvlapLen     [1/-1] bps. 
#     -rmOvlapRatio1  [-1] 0-1 (1 is 100%)
#     -rmOvlapRatio2  [-1] 0-1 (1 is 100%)
#     -rmOvlapType    [exon,CDS,mRNA]
#     -rmOvlapStrand  ['Both'] Could be 'Single'. 
#----------------------------------------------------------------------------------------------------
# -ovlapLongest     [Boolean] Keep only the longest overlapped genes. 
#                     Please note that I check only CDS' overlap if CDS feature exists. 
#                     If no CDS feature found, I will check for exon. 
#                     If no exon feature found, I will check for match_part. 
#                     And I will use the sum(length_of_features) instead of genomic_span as the topID's length. 
#   -ovlapFeatType    [CDS,exon,match_part] The order of features to be checked. Use the first existing one. 
#   -ovlapLength      [-1/1] The minimum length to say two genes are overlapped. 
#   -ovlapRatio       [-1] The minimum ratio to say two genes are overlapped. 
#                      Both -ovlapLength and -ovlapRatio will be checked when they are not '-1'; 
#                      If -ovlapRatio is '-1', -ovlapLength will be "1" by default. 
#   -ovlapStrand      ['Both'] Could be 'Single'
#----------------------------------------------------------------------------------------------------
# -islandGene       [4000] Given this flanking size, I provide genes with no other models in both flanking regions. 
#   -islandFeatType   [CDS,exon,match_part] The order of features to be checked. Use the first existing one. 
#                       Similar to -ovlapFeatType . 
#   -islandStrand     ['Both'] Could be 'Single'. 
#----------------------------------------------------------------------------------------------------
#
# 
# Filter options: 
# -extractFeat      [] Could be CDS. 
# -gff_top_hier     [undef()] Will use gffSunhh default value if not given. 
#                     Could be 'mrna,match,protein_match,expressed_sequence_match'
#
################################################################################
HH
	exit 1; 
}#sub usage() 


$opts{'help'} and &usage(); 
# Setting basic parameters: 
# Required options: -inGff
my $gff_obj = gffSunhh->new(); 
my $fas_obj = fastaSunhh->new(); 
my $mat_obj = mathSunhh->new(); 
my $iFh; 
if ( defined $opts{'inGff'} ) {
	$iFh = &openFH($opts{'inGff'}, '<'); 
} elsif ( @ARGV > 0 ) {
	$iFh = &openFH(shift, '<'); 
} elsif ( !( -t ) ) {
	$iFh = \*STDIN; 
} else {
	&usage(); 
}
my $oFh = \*STDOUT; 
defined $opts{'out'} and $oFh = &openFH($opts{'out'}, '>'); 
$opts{'seqInGff'} = ( $opts{'seqInGff'} ) ? 1 : 0; 
if (defined $opts{'gff_top_hier'}) {
	my @ta = grep { $_ ne '' } map { s!^\s+|\s+$!!g; lc($_); } split(/,/, $opts{'gff_top_hier'}); 
	$opts{'gff_top_hier'} = {}; 
	$gff_obj->_setTypeHash( $opts{'gff_top_hier'}, \@ta ); 
} else {
	$opts{'gff_top_hier'} = undef(); 
}
if ( defined $opts{'gffret'} ) {
	$opts{'idType'} //= 'mRNA'; 
}
if ( defined $opts{'idType'} ) { 
	$opts{'gff_top_hier'} = {}; 
	my @ta = grep { $_ ne '' } map { s!^\s+|\s+$!!g; lc($_); } split(/,/, $opts{'idType'}); 
	$gff_obj->_setTypeHash( $opts{'gff_top_hier'}, \@ta ); 
}
$opts{'sortGffBy'} //= 'lineNum'; 


# Step1. Read in gff file and sequence files. 
my (%in_gff, %in_seq ); 
{
	my ( $in_gff_href, $in_seq_href ) = $gff_obj->read_gff3File('gffFH'=>$iFh, 'saveFa'=>$opts{'seqInGff'}, 'top_hier'=>$opts{'gff_top_hier'} ); 
	%in_gff = %$in_gff_href; 
	%in_seq = %$in_seq_href; 
	if ( defined $opts{'scfFa'} ) {
		my $tmp_kseq_href = $fas_obj->save_seq_to_hash( 'faFile'=>$opts{'scfFa'} ); 
		for (keys %$tmp_kseq_href) {
			if (defined $in_seq{$_}) {
				&tsmsg("[Wrn] Skip repeated seqID [$_] in file [-scfFa $opts{'scfFa'}]\n"); 
			} else {
				( $in_seq{$_} = $tmp_kseq_href->{$_}{'seq'} ) =~ s!\s!!gs; 
			}
		}#End for 
	}
}

# Step2. Different actions. 
if ( $opts{'addFaToGff'} ) { 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'seq_href'=>\%in_seq, 'sort_by'=>$opts{'sortGffBy'} ); 
} elsif ( defined $opts{'compare2gffC'} ) {
	my $cFh = &openFH($opts{'compare2gffC'}, '<'); 
	my ($co_gff_href) = $gff_obj->read_gff3File('gffFH'=>$cFh, 'saveFa'=>0, 'top_hier'=>$opts{'gff_top_hier'}); 
	my %co_gff = %$co_gff_href; 
	undef($co_gff_href); 
	if ( defined $opts{'sameIntron'} ) {
		my ($kept_topIDs_aref) = &getID_sameIntron_gff(\%in_gff, \%co_gff, 'sameSingleExon'=>$opts{'sameSingleExon'}); 
		$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
	} elsif ( defined $opts{'rmOvlap'} ) {
		my ($kept_topIDs_aref) = &getID_rmOvlap_gff(
		 \%in_gff, \%co_gff, 
		 'rmOvlapLen'    => $opts{'rmOvlapLen'}, 
		 'rmOvlapRatio1' => $opts{'rmOvlapRatio1'}, 
		 'rmOvlapRatio2' => $opts{'rmOvlapRatio2'}, 
		 'rmOvlapType'   => $opts{'rmOvlapType'}, 
		 'rmOvlapStrand' => $opts{'rmOvlapStrand'}
		); 
		$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
	} else {
		&stopErr("[Err] Nothing to do with -compare2gffC\n"); 
	}
} elsif ( $opts{'ovlapLongest'} ) {
	$opts{'ovlapRatio'} //= -1; 
	$opts{'ovlapLength'} //= ( ( $opts{'ovlapRatio'} == -1 ) ? 1 : -1 ); 
	$opts{'ovlapStrand'} //= 'Both'; 
	$opts{'ovlapFeatType'} //= 'CDS,exon,match_part'; 
	my ($kept_topIDs_aref) = &getID_ovlapLongest_gff(
	 \%in_gff, 'ovlapLength'=>$opts{'ovlapLength'}, 
	 'ovlapRatio'=>$opts{'ovlapRatio'}, 
	 'ovlapStrand'=>$opts{'ovlapStrand'}, 
	 'ovlapFeatType'=>$opts{'ovlapFeatType'}
	); 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
} elsif ( defined $opts{'islandGene'} ) {
	$opts{'islandFeatType'} //= 'CDS,exon,match_part'; 
	$opts{'islandStrand'} //= 'Both'; 
	
	my ($kept_topIDs_aref) = &getID_islandGene_gff(
	 \%in_gff, 
	 'islandFeatType' => $opts{'islandFeatType'}, 
	 'islandStrand'   => $opts{'islandStrand'}, 
	 'islandFlank' => $opts{'islandGene'}
	); 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
} elsif ( defined $opts{'gffret'} ) {
	my $kept_topIDs_aref = []; 
	my $lisFh = &openFH($opts{'gffret'}, '<'); 
	while (<$lisFh>) {
		chomp; m!^\s*($|#)! and next; 
		my @ta = split(/\t/, $_); 
		push(@$kept_topIDs_aref, $ta[0]); 
	}
	close ($lisFh); 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
} elsif ( $opts{'listTopID'} ) {
	for ( sort { $in_gff{'ID2lineN'}{$a}<=>$in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
		print {$oFh} "$_\n"; 
	}
} else {
	&tsmsg("[Err] No valid action.\n"); 
	exit; 
}


# StepX. sub-routines here. 

# Input : (\%in_gff, '')
# Output: (\@topIDs_accepted)
sub getID_islandGene_gff {
	my $g1 = shift; 
	my %parm = $mat_obj->_setHashFromArr(@_); 
	$parm{'geneFeatType'} //= 'CDS,exon,match_part'; 
	$parm{'islandStrand'} //= 'Both'; 
	$parm{'islandStrand'} =~ m!^Both|Single$! or &stopErr("[Err] -islandStrand must be 'Both'/'Single'\n"); 
	$parm{'islandFlank'} //= 4000; 
	unless ($opts{'silent'}) {
		&tsmsg("[Rec][getID_islandGene_gff] Set -$_ as $parm{$_}\n") for ( sort keys %parm ); 
	}
	
	# Basic data
	my @back_topIDs; 
	my @chkTypes = split(/,/, $parm{'geneFeatType'}); 
	# Get locations. 
	my %chr2_topID_se_str_len; 
	my @topIDs = sort { $g1->{'ID2lineN'}{$a} <=> $g1->{'ID2lineN'}{$b} } keys %{ $g1->{'lineN_group'} }; 
	for my $topID ( @topIDs ) {
		my $chrID = &id_to_chr($g1, $topID); 
		my @arr; 
		for my $tmp_type ( @chkTypes ) {
			my %typeLocLis = &id_to_locLis($g1, $topID, [$tmp_type]); 
			if (@{$typeLocLis{'locLis'}} > 0) {
				@arr = ( $topID, $typeLocLis{'min'}, $typeLocLis{'max'}, $typeLocLis{'strand'}, $typeLocLis{'featLen'} ); 
				last; 
			}
		}# End for my $tmp_type (@chkTypes) 
		@arr > 0 or do { &tsmsg("[Wrn][getID_islandGene_gff] No countable locations found for $topID.\n"); next; }; 
		push(@{$chr2_topID_se_str_len{$chrID}}, [ @arr ]); 
	}
	
	if ( $parm{'islandStrand'} =~ m/^Both$/i ) {
		for my $chrID (keys %chr2_topID_se_str_len) {
			my ($ar) = &_rmClose_blk(
			 $chr2_topID_se_str_len{$chrID}, 
			 'flank'  => $parm{'islandFlank'}, 
			 'SEcol'  => [1,2]
			); 
			for my $tr (@$ar) {
				push(@back_topIDs, $tr->[0]); 
			}
		}#End for my $chrID
	} elsif ( $parm{'islandStrand'} =~ m/^Single$/i ) {
		for my $chrID (keys %chr2_topID_se_str_len) {
			my %tmp_str = map { $_->[3] => 1 } @{$chr2_topID_se_str_len{$chrID}}; 
			for my $t1 ( sort keys %tmp_str ) {
				my @t_ar = grep { $_->[3] eq $t1 } @{$chr2_topID_se_str_len{$chrID}}; 
				my ($t2) = &_rmClose_blk(
				 \@t_ar, 
				 'flank'  => $parm{'islandFlank'}, 
				 'SEcol'  => [1,2]
				); 
				for my $t3 (@$t2) {
					push(@back_topIDs, $t3->[0]); 
				}
			}
		}#End for my $chrID
	} else {
		&stopErr("[Err][getID_islandGene_gff] Unknown -islandStrand [$parm{'islandStrand'}]\n"); 
	}
	
	@back_topIDs = sort { $g1->{'ID2lineN'}{$a} <=> $g1->{'ID2lineN'}{$b} } @back_topIDs; 
	return (\@back_topIDs); 
}# sub getID_islandGene_gff () 


# Input : (\%in_gff, 'ovlapLength'=>[-1/1], 'ovlapRatio'=>[-1], 'ovlapStrand'=>'Both/Single', 'ovlapFeatType'=>'CDS,exon,match_part')
# Output: (\@topIDs_accepted)
sub getID_ovlapLongest_gff {
	my $g1 = shift; 
	my %parm = $mat_obj->_setHashFromArr(@_); 
	$parm{'ovlapRatio'} //= -1; 
	$parm{'ovlapLength'} //= ( ( $parm{'ovlapRatio'} == -1 ) ? 1 : -1 ); 
	$parm{'ovlapStrand'} //= 'Both'; 
	$parm{'ovlapStrand'} =~ m!^Both|Single$! or &stopErr("[Err] -ovlapStrand must be 'Both'/'Single'\n"); 
	$parm{'ovlapFeatType'} //= 'CDS,exon,match_part'; 
	unless ($opts{'silent'}) {
		&tsmsg("[Rec][getID_ovlapLongest_gff] Set -$_ as $parm{$_}\n") for (qw/ovlapRatio ovlapLength ovlapStrand ovlapFeatType/); 
	}
	
	# Basic data 
	my @back_topIDs; 
	my @ovlapTypes = split(/,/, $parm{'ovlapFeatType'}); 
	# Get locations. 
	my %chr2_topID_se_str_len; # {chrID} => [ [topID1, start1, end1, strand, featLen], [topID2, start2, end2, str, featLen], ... ]
	my @topIDs = sort { $g1->{'ID2lineN'}{$a} <=> $g1->{'ID2lineN'}{$b} } keys %{ $g1->{'lineN_group'} }; 
	for my $topID ( @topIDs ) {
		my $chr = &id_to_chr($g1, $topID); 
		my @arr; 
		for my $tmp_type (@ovlapTypes) {
			my %typeLocLis = &id_to_locLis($g1, $topID, [$tmp_type]); 
			if (@{$typeLocLis{'locLis'}} > 0) {
				#my @se = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } map { [ (sort { $a <=> $b } @{$_}[0,1]) , $_->[2] ] } @{$typeLocLis{'locLis'}}; 
				#my $len = 0; 
				#$len += ($_->[1]-$_->[0]+1) for (@se); 
				#@arr = ( $topID, $se[0][0], $se[-1][1], $se[0][2], $len ); 
				@arr = ( $topID, $typeLocLis{'min'}, $typeLocLis{'max'}, $typeLocLis{'strand'}, $typeLocLis{'featLen'} ); 
				last; 
			}
		}# End for my $tmp_type (@ovlapTypes) 
		@arr > 0 or do { &tsmsg("[Wrn][getID_ovlapLongest_gff] No countable locations found for $topID.\n"); next; }; 
		push(@{$chr2_topID_se_str_len{$chr}}, [ @arr ]); 
	}
	
	if ( $parm{'ovlapStrand'} =~ m/^Both$/i ) {
		for my $chrID (keys %chr2_topID_se_str_len) {
			# @{$chr2_topID_se_str_len{$chrID}} = sort { $b->[4] <=> $a->[4] || $a->[1] <=> $b->[1] || $a->[2]<=>$b->[2] } @{$chr2_topID_se_str_len{$chrID}}; 
			my ($ar) = &_rmOvlap_blk( $chr2_topID_se_str_len{$chrID}, 
			 'ovlLen'     => $parm{'ovlapLength'},  
			 'ovlRatio'   => $parm{'ovlapRatio'}, 
			 'ovlLenCol'  => 4, 
			 'ovlSEcol'   => [1,2], 
			 'ovlSrtCol'  => [4,1,2], 
			 'ovlSrtRule' => [-1,1,1], 
			); 
			for my $tr (@$ar) {
				push(@back_topIDs, $tr->[0]); 
			}
		}
	} elsif ( $parm{'ovlapStrand'} =~ m/^Single$/i ) {
		for my $chrID (keys %chr2_topID_se_str_len) {
			my %tmp_str = map { $_->[3] => 1 } @{$chr2_topID_se_str_len{$chrID}}; 
			for my $t1 (sort keys %tmp_str) {
				my @t_ar = grep { $_->[3] eq $t1 } @{$chr2_topID_se_str_len{$chrID}}; 
				my ($t2) = &_rmOvlap_blk( \@t_ar, 
				 'ovlLen'     => $parm{'ovlapLength'},  
				 'ovlRatio'   => $parm{'ovlapRatio'}, 
				 'ovlLenCol'  => 4, 
				 'ovlSEcol'   => [1,2], 
				 'ovlSrtCol'  => [4,1,2], 
				 'ovlSrtRule' => [-1,1,1], 
				); 
				for my $t3 (@$t2) {
					push(@back_topIDs, $t3->[0]); 
				}
			}
		}# End for my $chrID 
	} else {
		&stopErr("[Err] Unknown -ovlapStrand parameter.\n"); 
	}
	@back_topIDs = sort { $g1->{'ID2lineN'}{$a} <=> $g1->{'ID2lineN'}{$b} } @back_topIDs; 
	
	return (\@back_topIDs); 
}# sub getID_ovlapLongest_gff() 



# Input : (\@loc_arrays, 'flank'=>'4000', 'SEcol'=>[$col_S, $col_E]); 
# Output: (\@loc_arrays); 
sub _rmClose_blk {
	my $ar = shift; 
	my %parm = $mat_obj->_setHashFromArr(@_); 
	$parm{'flank'} //= 4000; 
	$parm{'SEcol'} //= [0,1]; 
	# Sort blocks. 
	my @srt = sort {
	 $a->[ $parm{'SEcol'}[0] ] <=> $b->[ $parm{'SEcol'}[0] ] || $a->[ $parm{'SEcol'}[1] ] <=> $b->[ $parm{'SEcol'}[1] ] 
	} @$ar; 
	# Detect good blocks. 
	my @back; 
	my %is_rm; 
	for ( my $i=0; $i<@srt; $i++ ) {
		my $is_good = 1; 
		my ( $cur_s, $cur_e ) = @{$srt[$i]}[ @{$parm{'SEcol'}} ]; 
		for ( my $j=$i+1; $j<@srt; $j++ ) {
			my ( $chk_s, $chk_e ) = @{$srt[$j]}[ @{$parm{'SEcol'}} ]; 
			$chk_s > ($cur_e+$parm{'flank'}) and last; 
			if ( $chk_s <= $cur_e+$parm{'flank'} ) {
				$is_rm{$i} = 1; 
				$is_rm{$j} = 1; 
				$is_good = 0; 
			}
		}# End for ( my $j=$i+1;
		$is_good == 1 and !(defined $is_rm{$i}) and push(@back, $srt[$i]); 
	}# for ( my $i=0; 
	
	return (\@back); 
}# _rmClose_blk



# Input  : (\@loc_arrays, 'ovlLen'=>'-1/1', 'ovlRatio'=>-1/0-1, 'ovlLenCol'=>ColN, 'ovlSEcol'=>[$col_S, $col_E], 'ovlSrtCol'=>[$col_S, $col_E], 'ovlSrtRule'=>[1, 1])
# Output : (\@loc_arrays); 
sub _rmOvlap_blk {
	my $ar = shift; 
	my %parm = $mat_obj->_setHashFromArr(@_); 
	$parm{'ovlRatio'}  //= -1; 
	$parm{'ovlRatio1'} //= $parm{'ovlRatio'} // -1; 
	$parm{'ovlRatio2'} //= $parm{'ovlRatio'} // -1; 
	unless (defined $parm{'ovlLen'}) {
		my $has_ratio = 0; 
		for ( qw/ovlRatio ovlRatio1 ovlRatio2/ ) {
			$parm{$_} == -1 and next; 
			$has_ratio = 1; 
			last; 
		}
		if ( $has_ratio == 1 ) {
			$parm{'ovlLen'} //= -1; 
		} else {
			$parm{'ovlLen'} //= 1; 
		}
	}
	$parm{'ovlSEcol'} //= [ 0,1 ]; 
	$parm{'ovlSrtCol'} //= $parm{'ovlSEcol'} ; 
	unless ( defined $parm{'ovlSrtRule'} ) {
		$parm{'ovlSrtRule'} = []; 
		for ( @{$parm{'ovlSrtCol'}} ) {
			push( @{$parm{'ovlSrtRule'}}, 1 ); 
		}
	}
	@{$parm{'ovlSrtCol'}} == @{$parm{'ovlSrtRule'}} or &stopErr("[Err][_rmOvlap_blk] ovlStrCol length is different from ovlSrtRule\n"); 
	# Sort array. 
	my @srt = sort {
	 my $return = 0; 
	 for ( my $i=0; $i<@{$parm{'ovlSrtCol'}}; $i++ ) {
	   $return  = $a->[ $parm{'ovlSrtCol'}[$i] ] <=> $b->[ $parm{'ovlSrtCol'}[$i] ] ; 
	   $return *= $parm{'ovlSrtRule'}[$i] ; 
	   $return and last; 
	 }
	 $return; 
	} @$ar; 
	# Record good blocks. 
	my @back; 
	my %is_rm; 
	for ( my $i=0; $i<@srt; $i++ ) {
		$is_rm{$i} and next; 
		push(@back, $srt[$i]); 
		my ( $cur_s, $cur_e ) = @{$srt[$i]}[ @{$parm{'ovlSEcol'}} ]; 
		$cur_s > $cur_e and ($cur_e, $cur_s) = ($cur_s, $cur_e); 
		my $cur_len = ( defined $parm{'ovlLenCol'} ) ? ($srt[$i][ $parm{'ovlLenCol'} ]) : ($cur_e-$cur_s+1) ; 
		for ( my $j=$i+1; $j<@srt; $j++ ) {
			$is_rm{$j} and next; 
			my ( $chk_s, $chk_e ) = @{$srt[$j]}[ @{$parm{'ovlSEcol'}} ]; 
			$chk_s > $chk_e and ($chk_e, $chk_s) = ($chk_s, $chk_e); 
			my $chk_len = ( defined $parm{'ovlLenCol'} ) ? ($srt[$j][ $parm{'ovlLenCol'} ]) : ($chk_e-$chk_s+1) ; 
			my ( $ovl_bp, $ovl_se ) = $mat_obj->ovl_region( $cur_s, $cur_e, $chk_s, $chk_e ); 
			$ovl_bp == 0 and next; 
			if ( $parm{'ovlLen'} > 0 ) {
				$ovl_bp >= $parm{'ovlLen'} and do { $is_rm{$j}=1; next; }; 
			}
			$parm{'ovlRatio1'} > 0 and $ovl_bp/$chk_len >= $parm{'ovlRatio1'} and do { $is_rm{$j}=1; next; }; 
			$parm{'ovlRatio2'} > 0 and $ovl_bp/$cur_len >= $parm{'ovlRatio2'} and do { $is_rm{$j}=1; next; }; 
		}# End for ( my $j=$i+1; 
	}# End for ( my $i=0; 
	
	return (\@back); 
}# _rmOvlap_blk


# Input : (\%in_gff_1, \%in_gff_2_Compare, 'rmOvlapLen'=>-1, 'rmOvlapRatio1'=>-1, 'rmOvlapRatio2'=>-1, 'rmOvlapType'=>'exon,CDS,match_part', 'rmOvlapStrand'=>'Both')
# Output: (\@back_topIDs); 
sub getID_rmOvlap_gff {
	my $g1 = shift; 
	my $g2 = shift; 
	my %parm = $mat_obj->_setHashFromArr(@_); 
	$parm{'rmOvlapRatio1'} //= -1; 
	$parm{'rmOvlapRatio2'} //= -1; 
	unless ( defined $parm{'rmOvlapLen'} ) {
		my $has_ratio = 0; 
		for ( qw/rmOvlapRatio1 rmOvlapRatio2/ ) {
			$parm{$_} == -1 and next; 
			$has_ratio = 1; 
			last; 
		}
		if ( $has_ratio == 1 ) {
			$parm{'rmOvlapLen'} //= -1; 
		} else {
			$parm{'rmOvlapLen'} //= 1; 
		}
	}#End unless 
	$parm{'rmOvlapType'} //= 'exon,CDS,match_part'; 
	$parm{'rmOvlapStrand'} //= 'Both'; 
	$parm{'rmOvlapStrand'} =~ m/^(Both|Single)$/i or &stopErr("[Err][getID_rmOvlap_gff] -rmOvlapStrand should be Both/Single.\n"); 
	
	my @back_topIDs; 
	my @chkTypes = split(/,/, $parm{'rmOvlapType'}); 
	
	my %chr2_topID_se_str_len_1; 
	my %chr2_topID_se_str_len_2; 
	my @topIDs_1 = sort { $g1->{'ID2lineN'}{$a} <=> $g1->{'ID2lineN'}{$b} } keys %{ $g1->{'lineN_group'} }; 
	my @topIDs_2 = sort { $g2->{'ID2lineN'}{$a} <=> $g2->{'ID2lineN'}{$b} } keys %{ $g2->{'lineN_group'} }; 
	for my $topID ( @topIDs_1 ) {
		my $chrID = &id_to_chr($g1, $topID); 
		my @arr; 
		for my $tmp_type ( @chkTypes ) {
			my %typeLocLis = &id_to_locLis($g1, $topID, [$tmp_type]); 
			if (@{$typeLocLis{'locLis'}} > 0) {
				@arr = ( $topID, $typeLocLis{'min'}, $typeLocLis{'max'}, $typeLocLis{'strand'}, $typeLocLis{'featLen'} ); 
				last; 
			}
		}# End for my $tmp_type (@chkTypes) 
		@arr > 0 or do { &tsmsg("[Wrn][getID_rmOvlap_gff] No countable locations found for $topID.\n"); next; }; 
		push(@{$chr2_topID_se_str_len_1{$chrID}}, [ @arr ]); 
	}
	for my $topID ( @topIDs_2 ) {
		my $chrID = &id_to_chr($g2, $topID); 
		my @arr; 
		for my $tmp_type ( @chkTypes ) {
			my %typeLocLis = &id_to_locLis($g2, $topID, [$tmp_type]); 
			if (@{$typeLocLis{'locLis'}} > 0) {
				@arr = ( $topID, $typeLocLis{'min'}, $typeLocLis{'max'}, $typeLocLis{'strand'}, $typeLocLis{'featLen'} ); 
				last; 
			}
		}# End for my $tmp_type (@chkTypes) 
		@arr > 0 or do { &tsmsg("[Wrn][getID_rmOvlap_gff] No countable locations found for $topID.\n"); next; }; 
		push(@{$chr2_topID_se_str_len_2{$chrID}}, [ @arr ]); 
	}
	
	if ( $parm{'rmOvlapStrand'} =~ m/^Both$/i ) {
		for my $chrID ( keys %chr2_topID_se_str_len_1 ) {
			if (!defined $chr2_topID_se_str_len_2{$chrID}) {
				&tsmsg("[Msg][getID_rmOvlap_gff] No chrID=$chrID found in co_gff\n"); 
				for my $tr_1 ( @{ $chr2_topID_se_str_len_1{$chrID} }  ) {
					push(@back_topIDs, $tr_1->[0]); 
				}
				next; 
			}
			my @ta_1 = sort { $a->[1]<=>$b->[1] || $a->[2]<=>$b->[2] } @{ $chr2_topID_se_str_len_1{$chrID} }; 
			my @ta_2 = sort { $a->[1]<=>$b->[1] || $a->[2]<=>$b->[2] } @{ $chr2_topID_se_str_len_2{$chrID} }; 
			for my $tr_1 ( @ta_1 ) {
				my $is_good = 1; 
				for my $tr_2 ( @ta_2 ) {
					$tr_2->[1] > $tr_1->[2] and last; 
					$tr_2->[2] < $tr_1->[1] and next; 
					my ($ar) = &_rmOvlap_blk( [ $tr_1, $tr_2 ], 
					 'ovlLen'     => $parm{'rmOvlapLen'},  
					 'ovlRatio1'  => $parm{'rmOvlapRatio1'}, 
					 'ovlRatio2'  => $parm{'rmOvlapRatio2'}, 
					 'ovlLenCol'  => 4, 
					 'ovlSEcol'   => [1,2], 
					 'ovlSrtCol'  => [1,2], 
					 'ovlSrtRule' => [1,1], 
					); 
					if (@$ar < 2) {
						$is_good = 0; 
						last; 
					}
				}
				$is_good == 1 and push( @back_topIDs, $tr_1->[0] ); 
			}
		}
	} elsif ( $parm{'rmOvlapStrand'} =~ m/^Single$/i ) {
		for my $chrID ( keys %chr2_topID_se_str_len_1 ) {
			if (!defined $chr2_topID_se_str_len_2{$chrID}) {
				&tsmsg("[Msg][getID_rmOvlap_gff] No chrID=$chrID found in co_gff\n"); 
				for my $tr_1 ( @{ $chr2_topID_se_str_len_1{$chrID} }  ) {
					push(@back_topIDs, $tr_1->[0]); 
				}
				next; 
			}
			
			my %tmp_str = map { $_->[3] => 1 } @{ $chr2_topID_se_str_len_1{$chrID} }; 
			for my $str (sort keys %tmp_str) {
				my @ta_1 = sort { $a->[1]<=>$b->[1] || $a->[2]<=>$b->[2] } grep { $_->[3] eq $str } @{ $chr2_topID_se_str_len_1{$chrID} }; 
				my @ta_2 = sort { $a->[1]<=>$b->[1] || $a->[2]<=>$b->[2] } grep { $_->[3] eq $str } @{ $chr2_topID_se_str_len_2{$chrID} }; 
				for my $tr_1 ( @ta_1 ) {
					my $is_good = 1; 
					for my $tr_2 ( @ta_2 ) {
						$tr_2->[1] > $tr_1->[2] and last; 
						$tr_2->[2] < $tr_1->[1] and next; 
						my ($ar) = &_rmOvlap_blk( [ $tr_1, $tr_2 ], 
						 'ovlLen'     => $parm{'rmOvlapLen'},  
						 'ovlRatio1'  => $parm{'rmOvlapRatio1'}, 
						 'ovlRatio2'  => $parm{'rmOvlapRatio2'}, 
						 'ovlLenCol'  => 4, 
						 'ovlSEcol'   => [1,2], 
						 'ovlSrtCol'  => [1,2], 
						 'ovlSrtRule' => [1,1], 
						); 
						if (@$ar < 2) {
							$is_good = 0; 
							last; 
						}
					}
					$is_good == 1 and push( @back_topIDs, $tr_1->[0] ); 
				}
			}
		}
	} else {
		&stopErr("[Err][getID_rmOvlap_gff] Unknown -rmOvlapStrand [$parm{'rmOvlapStrand'}]\n"); 
	}
	
	return (\@back_topIDs); 
}# sub getID_rmOvlap_gff() 

# Input : (\%in_gff_1, \%in_gff_2_Compare)
# Return: (\@topIDs_accepted)
sub getID_sameIntron_gff {
	my $g1 = shift; 
	my $g2 = shift; 
	my %parm = $mat_obj->_setHashFromArr(@_); 
	$parm{'sameSingleExon'} //= undef(); 
	
	my @topID_1 = sort { $g1->{'ID2lineN'}{$a} <=> $g1->{'ID2lineN'}{$b} } keys %{ $g1->{'lineN_group'} }; 
	my @topID_2 = sort { $g2->{'ID2lineN'}{$a} <=> $g2->{'ID2lineN'}{$b} } keys %{ $g2->{'lineN_group'} }; 
	my %chr2pid_2; 
	for my $topID (@topID_2) {
		my $chr = &id_to_chr($g2, $topID); 
		push(@{$chr2pid_2{$chr}}, $topID); 
	}
	# Output variables. 
	my @back_topIDs; 
	# Detect good topIDs. 
	my %intronLoc_infor_2; 
	CHK_TOPID_1: for my $tid1 ( @topID_1 ) {
		my ($chr1) = &id_to_chr($g1, $tid1); 
		my %exLocLis = &id_to_locLis($g1, $tid1, [qw/exon match_part/]); 
		my %intronLoc = &interval_locLis( $exLocLis{'locLis'}, 0 ); 
		CHK_TOPID_2: for my $tid2 ( @{$chr2pid_2{$chr1}} ) {
			my %intronLoc_2; 
			if (defined $intronLoc_infor_2{$tid2}) {
				%intronLoc_2 = %{ $intronLoc_infor_2{$tid2} }; 
			} else {
				my %exLocLis_2 = &id_to_locLis($g2, $tid2, [qw/exon match_part/]); 
				%intronLoc_2 = &interval_locLis( $exLocLis_2{'locLis'}, 0 ); 
				$intronLoc_infor_2{$tid2} = \%intronLoc_2; 
			}
			if ( $#{$intronLoc{'locLis'}} == -1 ) {
				# There is no intron in tid1; 
				if ( $#{$intronLoc_2{'locLis'}} == -1 ) {
					# tid2 has no intron either. 
					if ( $parm{'sameSingleExon'} ) {
						my ($is_same) = $mat_obj->compare_number_list( [ $intronLoc{'minmax'} ], [ $intronLoc_2{'minmax'} ], 'compare'=>'same' ); 
						if ( $is_same == 1 ) {
							push(@back_topIDs, $tid1); 
							next CHK_TOPID_1; 
						}
					} else {
						# Accept tid1 if tid1 is included by tid2
						my ( $spec1, $spec2 ) = $mat_obj->compare_number_list( [ $intronLoc{'minmax'} ], [ $intronLoc_2{'minmax'} ], 'compare'=>'nonovl' ); 
						if ( @$spec1 == 0 ) {
							push(@back_topIDs, $tid1); 
							next CHK_TOPID_1; 
						}
					}
				} else {
					# tid2 has intron, so not the same. 
					next CHK_TOPID_2; 
				}
			} else {
				# There are introns in tid1; 
				my ($is_same) = $mat_obj->compare_number_list( $intronLoc{'locLis'}, $intronLoc_2{'locLis'}, 'compare'=>'same' ); 
				if ( $is_same == 1 ) {
					push(@back_topIDs, $tid1); 
					next CHK_TOPID_1; 
				}
			}
		}#End for CHK_TOPID_2; 
	}#End for CHK_TOPID_1; 
	
	return (\@back_topIDs); 
}# sub getID_sameIntron_gff() 

# Input  : ([ [$s1, $e1], [$s2, $e2], ... ], with/without boundaries)
# Return : (%interval_LocLis)
#   {'locLis'} => [[$newS_1, $newE_1], [$newS_2, $newE_2], ...]
#   {'minmax'} => [$min, $max]
sub interval_locLis {
	my ($ar, $is_b) = @_; 
	$is_b //= 0; 
	my $add = ( $is_b ) ? 0 : 1 ; 
	my %back_lis; 
	$back_lis{'locLis'} = []; 
	my ($min, $max); 
	for my $tr (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] } @$ar) {
		my @se = sort { $a <=> $b } @{$tr}[0,1]; 
		$min //= $se[0]; 
		$max //= $se[1]; 
		$min > $se[0] and $min = $se[0]; 
		$max < $se[1] and $max = $se[1]; 
		if (@{$back_lis{'locLis'}} == 0) {
			push(@{$back_lis{'locLis'}}, [$se[1]+$add, ]); 
		} else {
			$back_lis{'locLis'}[-1][1] = $se[0]-$add; 
			push(@{$back_lis{'locLis'}}, [$se[1]+$add, ]); 
		}
	}
	pop(@{$back_lis{'locLis'}}); 
	$back_lis{'minmax'} = [$min, $max]; 
	return (%back_lis); 
}# sub interval_locLis() 

# Input  : (\%gff, $featID, [qw/exon match_part/])
# Return : (%locLis)
#   {scfID}   => scfID; 
#   {locLis}  => [ [exS_1, exE_1, strand(+/-)], [exS_2, exE_2, strand(+/-)], ... ]
#   {min}     => $minimum_of_exS
#   {max}     => $maximum_of_exE
#   {strand}  => '+'/'-', 
#   {featLen} => $number_length
sub id_to_locLis {
	my ( $gff, $featID, $accept_type ) = @_; 
	$accept_type //= ['exon', 'match_part']; 
	
	my $curLn = $gff->{'ID2lineN'}{$featID}; 
	my %back_lis = (); 
	$back_lis{'locLis'} = []; 
	my @off_lines_hash; 
	push( @off_lines_hash, map { $gff->{'lineN2hash'}{$_} } grep { defined $gff->{'lineN2hash'}{$_} } @{ $gff->{'lineN_group'}{$featID}{'curLn'} } ); 
	push( @off_lines_hash, map { $gff->{'lineN2hash'}{$_} } grep { defined $gff->{'lineN2hash'}{$_} } @{ $gff->{'lineN_group'}{$featID}{'offLn'} } ); 
	if ( @off_lines_hash == 0 ) {
		&tsmsg("[Wrn] No information found for featID [$featID]\n"); 
		return (%back_lis); 
	}
	($back_lis{'scfID'}) = &id_to_chr($gff, $featID); 
	for my $th ( @off_lines_hash ) {
		$th->{'seqID'} eq $back_lis{'scfID'} or do { &tsmsg("[Wrn] Skip line of differnt scfID [$th->{'seqID'}], which should be [$back_lis{'scfID'}].\n"); next; }; 
		my $is_good = 0; 
		for my $ctype ( @$accept_type ) {
			$th->{'type'} =~ m/^$ctype$/i and do { $is_good = 1; last; }; 
		}
		if ( $is_good == 1 ) {
			push(@{$back_lis{'locLis'}}, [ $th->{'start'}, $th->{'end'}, $th->{'strand'} ]); 
			$back_lis{'min'} //= $th->{'start'}; 
			$back_lis{'max'} //= $th->{'end'}; 
			$back_lis{'min'} > $th->{'start'} and $back_lis{'min'} = $th->{'start'}; 
			$back_lis{'max'} < $th->{'end'}   and $back_lis{'max'} = $th->{'end'}; 
			$back_lis{'strand'} //= $th->{'strand'}; 
			$back_lis{'featLen'} += ($th->{'end'}-$th->{'start'}+1); 
		}
	}
	return (%back_lis); 
}# sub id_to_locLis() 

# Input (\%in_gff, $featID)
# Return : ($chrID)
#   Usually it should be only one. 
#   Here I count for only the current and offspring features. 
sub id_to_chr {
	my ( $gff, $featID ) = @_; 
	my %chrID; 
	my $vv = 0; 
	# Check current feature. 
	my $curLn = $gff->{'ID2lineN'}{$featID}; 
	defined $gff->{'lineN2hash'}{$curLn} and $chrID{ $gff->{'lineN2hash'}{$curLn}{'seqID'} } = $vv; 
	## Check parents' feature. 
	#for my $pid ( keys %{$gff->{'PID2CID'}{$featID}} ) {
	#	my $parLn = $gff->{'ID2lineN'}{$pid}; 
	#	defined $gff->{'lineN2hash'}{$parLn} or next; 
	#	my $chr = $gff->{'lineN2hash'}{$parLn}{'seqID'}; 
	#	defined $chrID{$chr} and next; 
	#	$vv ++; 
	#	$chrID{$chr} = $vv; 
	#}
	# Check offsprings' feature. 
	for my $cid ( keys %{$gff->{'PID2CID'}{$featID}} ) {
		my $offLn = $gff->{'ID2lineN'}{$cid}; 
		defined $gff->{'lineN2hash'}{$offLn} or next; 
		my $chr = $gff->{'lineN2hash'}{$offLn}{'seqID'}; 
		defined $chrID{ $chr } and next; 
		$vv ++; $chrID{ $chr } = $vv; 
		last; 
	}
	(keys %chrID) > 0 or &stopErr("[Err] No chrID found for featID [$featID]\n"); 
	my @chrIDs = sort { $chrID{$a}<=>$chrID{$b} } keys %chrID; 
	if (@chrIDs > 1) {
		&tsmsg("[Err] Found multiple chrIDs for featID [$featID]: @chrIDs\n"); 
	}
	return ($chrIDs[0]); 
}#sub id_to_chr() 


