#!/usr/bin/perl
# 2016-07-08 Add function to extract CDS sequence according to gff file. 
# 2018-05-14 Add function to extract CDS sequence according to gff file. 
# 2018-07-10 Sort GFF by its inner features. 
# 2019-01-16 Change 'frame' value to fit blastx and transeq; 
# 2021-07-09 Fix 'frame' to fit deal_fasta.pl. 
# 2022-01-19 Add gff_top_hier definition in opts variable.
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
	"sortGffBy:s", # Could be 'raw/lineNum/str_posi/position'
	"silent!", 
	"outFas:s",    # Output fasta file if defined. 
	
	# Actions 
	"sort!", 
	"simpleSort!", 
	"seqret!", # Not added. 
	  "extractFeat:s", # Not used. 

	"getAGP:s", 

	"fixTgt!", # 
	
	"compare2gffC:s", 
	 "sameAll!", # Not added. 
	 "sameIntron!", "sameSingleExon!", # Finished. 
	 "rmOvlap!", "rmOvlapLen:i", "rmOvlapRatio1:f", "rmOvlapRatio2:f", "rmOvlapType:s", "rmOvlapStrand:s", # To be improved. 
	
	"ovlapLongest!", "ovlapLength:i", "ovlapRatio:f", "ovlapStrand:s", "ovlapFeatType:s", # Finished. 
	"islandGene:i", "islandFeatType:s", "islandStrand:s", # Finished. 
	"gffret:s", "idType:s", # FInished. 
	"gffInRegion:s", "bothStr!", "onlyFull!", # 
	"listTopID!", # Finished. 
	"getLoc:s", "joinLoc!", # Doing. Could be mRNA/CDS/gene/intron, others may work but not guaranteed. 
	"getJnLoc!", # Doing. 
	
	"addFaToGff!", # Finished. 
	
	"list_intron!",  # Finished. 
	 "intron_byFeat:s", # [CDS], could be 'CDS|Exon'. 
	
	"ch_makerID:s", # Finished. The input gff must be exactly as maker output. 
	 "geneID_list:s", # [file], "oldGeneID \t newGeneID"
	
	"ch_ID:s",      # Change feature IDs in GFF file; 

	"ch_locByAGP:s", # Finished. I think I should add 'sortTopIDBy' to &write_gff3File() . 
	 "sortTopIDBy:s", # raw/lineNum/position
	
	# Filter options. 
	"gff_top_hier:s", # Could be 'mrna,match,protein_match,expressed_sequence_match'
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
# -outFas           [] Not output by default. Not fully supported yet. 
# 
# Action options: 
#----------------------------------------------------------------------------------------------------
# -simpleSort       [Boolean]
#----------------------------------------------------------------------------------------------------
# -ch_makerID       [oldID_newID] Convert mRNA ID and gene IDs according to file 'oldID_newID'. 
#                     'oldID_newID' should have two column, old ID and new ID. 
#   -geneID_list    [oldID_newID] Additionally provide a list to change gene IDs. 
#----------------------------------------------------------------------------------------------------
# -ch_ID            [oldID_newID] All the ID/Parent should exist in the input oldID_newID list. 
#----------------------------------------------------------------------------------------------------
# -ch_locByAGP      [in_ref.agp] Change IDs and locations according to this agp. Here assume to change agp_ctg to agp_scaff. 
#   -sortTopIDBy    ['raw'] raw/lineNum/position. In fact, raw should be equal to lineNum. 
#----------------------------------------------------------------------------------------------------
# -addFaToGff       [Boolean] Add -scfFa to the tail of -inGff . 
#----------------------------------------------------------------------------------------------------
# -seqret           [Boolean] Retrieve fasta sequence. Not used yet. 
#   -extractFeat    [] Could be CDS. Not used yet. 
#
#----------------------------------------------------------------------------------------------------
# -fixTgt           [Boolean] Fix coordinates problem in "Target:" section. 
#                     Change "Target=ID 0 100" to "Target=ID 1 100"
#
#----------------------------------------------------------------------------------------------------
# -gffret           [input_ID_list] I use the first column. 
#   -idType         [mRNA] Could be 'gene/exon/CDS/match/match_part'
#----------------------------------------------------------------------------------------------------
# -gffInRegion      [input_Region_list] Format: SeqID \\t Start \\t End [\\t Strand(+/-/.)]
#                     If defined Strand column, only genes overlapping in the same strand will be extracted. 
#   -idType         ['mRNA'] Decide the feature type which is checked for overlapping. 
#   -bothStr        [Boolean] Ignore strand information if given. 
#   -onlyFull       [Boolean] Only fully included genes will be extracted. 
#----------------------------------------------------------------------------------------------------
# -listTopID        [Boolean] 
#   -idType         [Same to -gff_top_hier if not given]. 
#----------------------------------------------------------------------------------------------------
# -getLoc           ['mRNA'] Could be a feature of mRNA/CDS/gene/intron, others may work but not guaranteed. 
#                     Output a table in format : 
#                       feat_ID\\tfeatParent_ID\\tchr_ID\\tfeat_Start\\tfeat_End\\tStrand(+/-)
#                     In which "feat_ID" is defined in feature. 
#                     If no featParent_ID exists, set featParent_ID=featID
#   -joinLoc        [Boolean] Output joined table by featParent_ID, the format is like : 
#                     LengthInFeat\\tfeatParent_ID\\tchrID\\tfeat_Start_min\\tfeat_End_max\\tStrand(+/-)\tfeat_Start_1,feat_End_1;feat_Start2,feat_End_2;...
#----------------------------------------------------------------------------------------------------
# -getJnLoc         [Boolean]
#----------------------------------------------------------------------------------------------------
# -list_intron      [Boolean] Provide intron table according to '-intron_byFeat'. 
#   -intron_byFeat  [CDS] Could be a feature of CDS/exon ; 
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
#     -rmOvlapLen     [1/-1] bps. Value '-1' means not considering by overlapped length. 
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
# -getAGP           ['CDS'] Get .agp file for CDS features; Not used yet. 
#----------------------------------------------------------------------------------------------------
#
# 
# Filter options: 
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
my $oFasFh = undef(); 
defined $opts{'outFas'} and $oFasFh = &openFH($opts{'outFas'}, '>'); 

if ( $opts{'getJnLoc'} ) {
	&action_getJnLoc(); 
	exit(); 
} elsif ( defined $opts{'ch_ID'} ) {
	&action_ch_ID(); 
	exit(); 
} elsif ( defined $opts{'simpleSort'} ) {
	&action_simpleSort(); 
	exit(); 
}

sub action_simpleSort {
	my %trans_feature = qw(
		contig             contig
		gene               gene
		transcript         mrna
		mrna               mrna
		five_prime_utr     five_prime_utr
		exon               exon
		cds                cds
		three_prime_utr    three_prime_utr
	); 
	my %good_feature = qw(
		contig             0.5
		gene               1
		mrna               2
		five_prime_utr     3
		exon               4 
		cds                5 
		three_prime_utr    6
	); 
	my (@lines, %info, $lineN); 
	$lineN = 0; 
	while (<$iFh>) {
		$lineN ++; 
		chomp; 
		if (m!^\s*#!) {
			if (@lines > 0) {
				# Sort and output @lines; 
				@lines = (); 
			}
			print {$oFh} "$_\n"; 
			next; 
		}
		my @ta = split(/\t/, $_); 
		my $featID = lc($ta[2]); 
		defined $trans_feature{$featID} or &stopErr("[Err] Unknown feature [$featID]\n"); 
		$info{'ln2line'}{$lineN} = $_; 
		$info{'ln2feat'}{$lineN} = $trans_feature{$featID}; 
		if ($ta[8] =~ m!(?:^|;)\s*Parent\s*=\s*([^\s;]+)\s*(?:$|;)!i) {
			my $pID = $1; 
			$pID =~ m!,! and &stopErr("[Err] I don't support multiple parents.\n"); 
			$info{'ln2parentID'}{$lineN} = $pID; 
		}
		if ($ta[8] =~ m!(?:^|;)\s*ID\s*=\s*([^\s;]+)\s*(?:$|;)!i) {
			my $cID = $1; 
			$info{'ln2currentID'}{$lineN} = $cID; 
		}
	}
}# action_simpleSort() 



$opts{'seqInGff'} = ( $opts{'seqInGff'} ) ? 1 : 0; 
if (defined $opts{'gff_top_hier'}) {
	my @ta = grep { $_ ne '' } map { s!^\s+|\s+$!!g; lc($_); } split(/,/, $opts{'gff_top_hier'}); 
	$opts{'gff_top_hier'} = {}; 
	$gff_obj->_setTypeHash( $opts{'gff_top_hier'}, \@ta ); 
} else {
	$opts{'gff_top_hier'} = undef(); 
}
if ( defined $opts{'gffret'} or defined $opts{'gffInRegion'} ) {
	$opts{'idType'} //= 'mRNA'; 
}
if ( defined $opts{'idType'} ) { 
	$opts{'gff_top_hier'} = {}; 
	my @ta = grep { $_ ne '' } map { s!^\s+|\s+$!!g; lc($_); } split(/,/, $opts{'idType'}); 
	$gff_obj->_setTypeHash( $opts{'gff_top_hier'}, \@ta ); 
}
if ( defined $opts{'getLoc'} ) {
	$opts{'gff_top_hier'} = {}; 
	my @ta = grep { $_ ne '' } map { s!^\s+|\s+$!!g; lc($_); } split(/,/, $opts{'getLoc'}); 
	$gff_obj->_setTypeHash( $opts{'gff_top_hier'}, \@ta ); 
}
if ( defined $opts{'list_intron'} ) {
	$opts{'gff_top_hier'} = {}; 
	$opts{'intron_byFeat'} //= 'CDS'; 
	$gff_obj->_setTypeHash( $opts{'gff_top_hier'}, [$opts{'intron_byFeat'}] ); 
}
if ( defined $opts{'getAGP'} ) {
	if ( $opts{'getAGP'} =~ m!^CDS$!i ) {
		$opts{'gff_top_hier'} = { 'mrna'=>1, 'match'=>2, 'protein_match'=>3, 'expressed_sequence_match'=>4 }; 
		$opts{'seqret'} = 1; 
		$opts{'extractFeat'} = 'CDS'; 
	} elsif ( $opts{'getAGP'} =~ m!^exon$!i ) {
		$opts{'gff_top_hier'} = { 'mrna'=>1, 'match'=>2, 'protein_match'=>3, 'expressed_sequence_match'=>4 }; 
		$opts{'seqret'} = 1; 
		$opts{'extractFeat'} = 'exon'; 
	} else {
		&stopErr("[Err] Unknown type for -getAGP ; Should be CDS/exon\n"); 
	}
}
if ( $opts{'seqret'} ) {
	if ( $opts{'extractFeat'} =~ m!^(cds|exon)$!i) {
		$opts{'gff_top_hier'} = { 'mrna'=>1, 'match'=>2, 'protein_match'=>3, 'expressed_sequence_match'=>4 }; 
	}
}

$opts{'sortGffBy'} //= 'lineNum'; 
$opts{'sortTopIDBy'} //= 'raw'; $opts{'sortTopIDBy'} = lc($opts{'sortTopIDBy'}); # raw/lineNum/position 
$opts{'extractFeat'} //= 'CDS'; 

# Step1. Read in gff file and sequence files. 
my (%in_gff, %in_seq ); 
&load_gff_fas( \%in_gff, \%in_seq ); 

# Step2. Different actions. 
if ( $opts{'sort'} ) {
	&action_sort(); 
} elsif ( $opts{'addFaToGff'} ) { 
	&action_addFaToGff(); 
} elsif ( defined $opts{'compare2gffC'} ) {
	&action_compare2gffC(); 
} elsif ( $opts{'gffInRegion'} ) {
	&action_gffInRegion(); 
} elsif ( $opts{'ovlapLongest'} ) {
	&action_ovlapLongest(); 
} elsif ( defined $opts{'islandGene'} ) {
	&action_islandGene(); 
} elsif ( defined $opts{'gffret'} ) {
	&action_gffret(); 
} elsif ( $opts{'listTopID'} ) {
	&action_listTopID(); 
} elsif ( defined $opts{'getLoc'} ) {
	&action_getLoc(); 
} elsif ( $opts{'fixTgt'} ) { 
	&action_fixTgt(); 
} elsif ( $opts{'list_intron'} ) {
	&action_list_intron(); 
} elsif ( defined $opts{'ch_makerID'} ) {
	&action_ch_makerID(); 
} elsif ( defined $opts{'ch_locByAGP'} ) {
	&action_ch_locByAGP(); 
} elsif ( $opts{'seqret'} ) {
	&action_seqret(); 
} else {
	&tsmsg("[Err] No valid action.\n"); 
	exit; 
}

################################################################################
############### StepX. action sub-routines here. 
################################################################################
sub action_seqret {
	my %str2num = qw(
	 +     1
	 -     -1
	 1     1
	 -1    -1
	 0     -1
	 plus  1
	 minus -1
	); 
	TOPID:
	for my $topID ( sort { $in_gff{'ID2lineN'}{$a} <=> $in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
		my @posi_top; 
		my @posi_cds; 
		my $curLnNum = $in_gff{'lineN_group'}{$topID}{'curLn'}[0]; 
		my $top_str = ''; 
		my $top_chr = ''; 
		my $top_name = ''; 
		my %curLnH; 
		if ( defined $in_gff{'lineN2hash'}{$curLnNum} ) {
			%curLnH = %{ $in_gff{'lineN2hash'}{$curLnNum} }; 
			$top_str = $str2num{ lc( $curLnH{'strand'} ) } // ''; 
			$top_chr = $curLnH{'seqID'} // ''; 
			$top_name = $curLnH{'attrib'}{'featID'} // ''; 
			@posi_top = ( [ $curLnH{'start'}, $curLnH{'end'} ] ); 
		}
		OFFID:
		for my $offLnNum ( @{$in_gff{'lineN_group'}{$topID}{'offLn'}} ) {
			defined $in_gff{'lineN2hash'}{$offLnNum} or next OFFID; 
			my %offLnH = %{ $in_gff{'lineN2hash'}{$offLnNum} }; 
			if ( $opts{'extractFeat'} =~ m!^cds$!i ) {
				$offLnH{'type'} =~ m!^CDS$!i or next OFFID; 
			} elsif ( $opts{'extractFeat'} =~ m!^exon$!i ) {
				$offLnH{'type'} =~ m!^exon$!i or next OFFID; 
			}
			$top_str eq '' and $top_str = $str2num{ lc($offLnH{'strand'}) } // ''; 
			$top_chr eq '' and $top_chr = $offLnH{'seqID'} // ''; 
			if ( keys %{$offLnH{'attrib'}{'parentID'}} > 0 ) {
				my @ta1 = sort keys %{ $offLnH{'attrib'}{'parentID'} }; 
				my $ta1_txt = join(";", @ta1); 
				if ( $top_name eq '' ) {
					$top_name = $ta1_txt; 
				} else {
					$top_name eq $ta1_txt or &stopErr("[Err] Unequal top_name : [ $top_name, $ta1_txt ]\n"); 
				}
			}
			my @ta = split(/\t/, $in_gff{'lineN2line'}{$offLnNum}); 
			my $t_frame; 
			if ( $ta[7] eq '.' or $ta[7] eq '0' ) {
				$t_frame = 1; 
			} elsif ( $ta[7] == 1 ) {
				$t_frame = 3; 
			} elsif ( $ta[7] == 2 ) {
				$t_frame = 2; 
			} else {
				&stopErr("[Err] Failed to parse gff phase [$ta[7]]\n"); 
			}
			push( @posi_cds, [ $offLnH{'start'}, $offLnH{'end'}, $t_frame] ); 
		}# End for my $offLnNum 
		@posi_cds > 0 or next TOPID; 
		$top_str eq '' and do { &tsmsg("[Wrn] No strand information for topID=[$topID]\n"); $top_str = 1; }; 
		$top_chr eq '' and do { &stopErr("[Err] No top_chr found for [$topID]\n"); }; 
		defined $in_seq{$top_chr} or do { &tsmsg("[Wrn] No [$top_chr] sequence found for [$topID]\n"); next TOPID; }; 
		# Check @posi_cds
		for my $tr (@posi_cds) {
			$tr->[0] > $tr->[1] and &stopErr("[Err] I can't accept start[$tr->[0]] > end[$tr->[1]] in gff3 file.\n"); 
		}
		if ($top_str == -1) {
			# Be careful when the CDS loci should not ordered by position; 
			@posi_cds = sort { $b->[0] <=> $a->[0] } @posi_cds; 
		}
		if ( defined $opts{'getAGP'} ) {
			my %h1; 
			$h1{'eleN'} = 0; 
			# ##agp-version   2.0
			# Super_scaffold_248      1       341379  1       W       SpoScf_00639    1       341379  -
			# Super_scaffold_248      341380  600696  2       N       259317  scaffold        yes     map
			for my $tr (@posi_cds) {
				$h1{'eleN'} ++; 
				$h1{'prevE'} //= 0; 
				print {$oFh} join("\t", 
					$top_name, 
					$h1{'prevE'}+1, 
					$h1{'prevE'}+$tr->[1]-$tr->[0]+1, 
					$h1{'eleN'}, 
					'W', 
					$top_chr, 
					$tr->[0], 
					$tr->[1], 
					( $top_str == -1 ) ? '-' : '+'
				)."\n"; 
				$h1{'prevE'} = $h1{'prevE'}+$tr->[1]-$tr->[0]+1; 
			}
			next TOPID; 
		}

		# get sequences; 
		my @sub_seqs; 
		for my $tr ( @posi_cds ) {
			push( @sub_seqs, substr($in_seq{$top_chr}, $tr->[0]-1, $tr->[1]-$tr->[0]+1) ); 
		}
		$top_str == -1 and @sub_seqs = &rev_comp(@sub_seqs); 
		my $final_seq = join('', @sub_seqs); 
		if ( $opts{'extractFeat'} =~ m!^cds$!i ) {
			print {$oFh} ">$top_name [frame=$posi_cds[0][2]]\n$final_seq\n"; 
		} elsif ( $opts{'extractFeat'} =~ m!^exon$!i ) {
			print {$oFh} ">$top_name\n$final_seq\n"; 
		}
	}# End for my $topID 
}# action_seqret() 

sub action_getJnLoc {
	my $wrk_dir = &fileSunhh::new_tmp_dir(); 
	mkdir($wrk_dir); 
	open O,'>',"$wrk_dir/in.gff" or die; 
	while (<$iFh>) {
		print O $_; 
	}
	close O; 
	open I1,'-|',"perl $0 -getLoc mRNA -inGff $wrk_dir/in.gff" or die; 
	open O1,'>',"$wrk_dir/in.gff.loc_mrna" or die; 
	while (&wantLineC( \*I1 )) {
		if ($. == 1) {
			my @ta = &splitL("\t", $_); 
			$ta[0] eq 'FeatID' and $ta[0] = 'mrnaID'; 
			$ta[1] eq 'ParentID' and $ta[1] = 'geneID'; 
			@ta[3,4,5] = qw/mrnaStart mrnaEnd mrnaStrand/; 
			print O1 join("\t", @ta)."\n"; 
			next; 
		}
		print O1 "$_\n"; 
	}
	close O1; 
	close I1; 
	open I2,'-|',"perl $0 -getLoc CDS -joinLoc -inGff $wrk_dir/in.gff" or die; 
	open O2,'>',"$wrk_dir/in.gff.loc_cds" or die; 
	while ( &wantLineC( \*I2 ) ) {
		if ($. == 1) {
			my @ta = &splitL("\t", $_); 
			@ta = qw/LenInCDS       mrnaID  SeqID   CDSStart     CDSEnd    CDSStrand       CDSBlocks/; 
			print O2 join("\t", @ta)."\n"; 
			next; 
		}
		print O2 "$_\n"; 
	}
	close I2; 
	close O2; 
	open I3,'-|', "ColLink.pl $wrk_dir/in.gff.loc_cds -f1 $wrk_dir/in.gff.loc_mrna -keyC1 0 -keyC2 1 -add -COl1 1,3,4,5 | deal_table.pl -column 1,7,2,8,9,10,3,4,0,6" or die; 
	while (<I3>) {
		chomp($_); 
		my @ta = &splitL("\t", $_); 
		if ($. == 1) {
			$ta[0] eq 'mrnaID' or &stopErr("[Err] Bad first line: $_\n"); 
			print {$oFh} join("\t", @ta, qw/CDSBlocksNum 5UTR 3UTR/)."\n"; 
			next; 
		}
		my $cdsN = ( $ta[9] =~ tr/;/;/ ) + 1; 
		my ($utr5, $utr3); 
		if ($ta[5] eq '-') {
			$utr5 = $ta[4]-$ta[7]; 
			$utr3 = $ta[6]-$ta[3]; 
		} elsif ($ta[5] =~ m!^(\.|\+)$!) {
			$utr5 = $ta[6]-$ta[3]; 
			$utr3 = $ta[4]-$ta[7]; 
		} else {
			&stopErr("[Err] Bad strand character [$ta[5]] in line:\n$_\n"); 
		}
		print {$oFh} join("\t", @ta, $cdsN, $utr5, $utr3)."\n"; 
	}
	close I3; 

	&fileSunhh::_rmtree($wrk_dir); 

	return; 

}# action_getJnLoc() 
sub action_ch_locByAGP {
	my %info_ctg2scf; # {ctgID}=>[ [ctgS, ctgE, scfID, scfS, scfE, scfStr(+/-/?)], [], ... ]
	%info_ctg2scf = %{ &fileSunhh::load_agpFile( $opts{'ch_locByAGP'} ) }; 
	
	for my $lineN (sort {$a <=> $b} keys %{$in_gff{'lineN2line'}}) {
		defined $in_gff{'lineN2hash'}{$lineN} or next; 
		my $line_txt = $in_gff{'lineN2line'}{$lineN}; 
		my @ta = &splitL("\t", $line_txt); 
		my ($ctgID, $ctgS, $ctgE, $ctgStr) = @ta[0, 3,4, 6]; 
		defined $info_ctg2scf{$ctgID} or next; 
		my @loc1 = $mat_obj->switch_position( 'qry2ref' => \%info_ctg2scf, 'qryID'=>$ctgID, 'qryPos'=>$ctgS, 'qryStr'=>$ctgStr ); 
		my @loc2 = $mat_obj->switch_position( 'qry2ref' => \%info_ctg2scf, 'qryID'=>$ctgID, 'qryPos'=>$ctgE, 'qryStr'=>$ctgStr ); 
		( @loc1 > 1 or @loc2 > 1 ) and do { &tsmsg("[Wrn] Multiple loci found for ($ctgID, $ctgS, $ctgE, $ctgStr)\n"); }; 
		$loc1[0][0] eq $loc2[0][0] or &stopErr("[Err] Found different scaffold IDs for ($ctgID, $ctgS, $ctgE, $ctgStr): $loc1[0][0] & $loc2[0][0]\n"); 
		$ta[0] = $loc1[0][0]; 
		$ta[3] = $loc1[0][1]; 
		$ta[4] = $loc2[0][1]; 
		$ta[6] = $loc1[0][2]; 
		@ta[3,4] = sort {$a<=>$b} @ta[3,4]; 
		
		$in_gff{'lineN2line'}{$lineN} = join("\t", @ta); 
		# Replace 'lineN2hash' 
		$in_gff{'lineN2hash'}{$lineN}{'seqID'}  = $ta[0]; 
		$in_gff{'lineN2hash'}{$lineN}{'start'}  = $ta[3]; 
		$in_gff{'lineN2hash'}{$lineN}{'end'}    = $ta[4]; 
		$in_gff{'lineN2hash'}{$lineN}{'strand'} = $ta[6]; 
	}
	
	# Give topID loc: 
	my %topID_loc; 
	my @topIDs; 
	my %doneID; 
	LINE_N: 
	for my $lineN (sort {$a <=> $b} keys %{$in_gff{'lineN2line'}}) {
		defined $in_gff{'lineN2hash'}{$lineN} or next; 
		my $lineN2hash_href = $in_gff{'lineN2hash'}{$lineN}; 
		
		my $featID = $in_gff{'lineN2ID'}{$lineN}; 
		defined $doneID{$featID} and next LINE_N; 
		defined $in_gff{'lineN_group'}{$featID} and do { $topID_loc{$featID} = [ @{$lineN2hash_href}{qw/seqID start end strand/} ]; next LINE_N;  }; 
		
		my @parentIDs = keys %{$in_gff{'CID2PID'}{$featID}}; 
		for my $pID (@parentIDs) {
			defined $doneID{$pID} and next; 
			defined $in_gff{'lineN_group'}{$pID} or next; 
			my $p_lineN = $in_gff{'ID2lineN'}{$pID}; 
			if (defined $in_gff{'lineN2hash'}{$p_lineN}) {
				$topID_loc{$pID} = [ @{$in_gff{'lineN2hash'}{$p_lineN}}{qw/seqID start end strand/} ]; 
			} else {
				$topID_loc{$pID} = [ @{$lineN2hash_href}{qw/seqID start end strand/} ]; 
			}
			$doneID{$pID} = 1; 
		}
	}
	
	
	my $use_topID = 0; 
	my $lineN_group_href = $in_gff{'lineN_group'}; 
	if ($opts{'sortTopIDBy'} eq 'raw') { 
		# # raw/lineNum/position 
		$use_topID = 0 ; 
	} elsif ($opts{'sortTopIDBy'} eq 'linenum') {
		$use_topID = 1 ; 
		@topIDs = sort { $lineN_group_href->{$a}{'curLn'}[0]<=>$lineN_group_href->{$b}{'curLn'}[0] } keys %$lineN_group_href ; 
	} elsif ($opts{'sortTopIDBy'} eq 'position') {
		$use_topID = 1 ; 
		for (keys %$lineN_group_href) {
			defined $topID_loc{$_} or die "$_\n"; 
		}
		@topIDs = sort { 
		 $topID_loc{$a}[0] cmp $topID_loc{$b}[0] 
		 || $topID_loc{$a}[1] <=> $topID_loc{$b}[1]
		 || $topID_loc{$a}[2] <=> $topID_loc{$b}[2]
		 || $topID_loc{$a}[3] cmp $topID_loc{$b}[3]
		} keys %$lineN_group_href; 
	} else {
		&tsmsg("[Wrn] Unknown -sortTopIDBy [$opts{'sortTopIDBy'}] skipped\n"); 
	}

	$use_topID == 0 and $gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'sort_by'=>$opts{'sortGffBy'} ); 
	$use_topID == 1 and $gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'sort_by'=>$opts{'sortGffBy'}, 'topIDs_aref' => \@topIDs ); 
	
	return(); 
}# sub action_ch_locByAGP() 
sub action_ch_ID {
	my %o2n = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $opts{'ch_ID'} ); 
	while (<$iFh>) {
		chomp; 
		if (m!^\s*(#|$)!) {
			print {$oFh} "$_\n"; 
			next; 
		}
		my @ta = split(/\t/, $_); 
		if ($ta[8] =~ m!ID=([^\s;]+)!) {
			my $id0=$1; 
			my @id1 = split(/,/, $id0); 
			for my $id2 (@id1) {
				defined $o2n{$id2} and $id2 = $o2n{$id2}; 
			}
			$id0=join(',', @id1);
			$ta[8] =~ s!ID=([^\s;]+)!ID=$id0!; 
		}
		if ($ta[8] =~ m!Parent=([^\s;]+)!) {
			my $id0=$1; 
			my @id1 = split(/,/, $id0); 
			for my $id2 (@id1) {
				defined $o2n{$id2} and $id2 = $o2n{$id2}; 
			}
			$id0=join(',', @id1);
			$ta[8] =~ s!Parent=([^\s;]+)!Parent=$id0!; 
		}
		print {$oFh} join("\t", @ta)."\n"; 
	}
}# action_ch_ID() 
sub action_ch_makerID {
	my $fh = &openFH( $opts{'ch_makerID'}, '<' ); 
	my %o2n_m; 
	my %o2n_g; 
	# my %n2o_g; 
	while (<$fh>) {
		&isSkipLine($_) and next; 
		chomp; s/[^\t\S]+$//; 
		my @ta = &splitL("\t", $_); 
		$o2n_m{$ta[0]} //= $ta[1]; 
		$ta[0] =~ s/\-mRNA\-\d+$//i; 
		$o2n_g{$ta[0]} //= "$ta[1]-gene"; 
		# push(@{$n2o_g{ "$ta[1]-gene" }}, $ta[0] ); 
	}
	close($fh); 
	if ( defined $opts{'geneID_list'} ) {
		my $fh_g = &openFH( $opts{'geneID_list'}, '<' ); 
		while (<$fh_g>) {
			&isSkipLine($_) and next; 
			chomp; s/[^\t\S]+$//; 
			my @ta = &splitL("\t", $_); 
			$o2n_g{$ta[0]} = $ta[1]; 
		}
		close($fh_g); 
	}
	
	for my $line_txt (values %{$in_gff{'lineN2line'}}) {
		my @ta = &splitL("\t", $line_txt); 
		scalar(@ta) >= 9 or next; 
		my %attrHash = %{ $gff_obj->_getAttrHash( 'attribText'=>$ta[8] ) }; 
		if ( lc($ta[2]) eq 'gene' ) {
			if ( defined $o2n_g{$attrHash{'featID'}} ) {
				my $s=$attrHash{'featID'}; 
				my $v=$o2n_g{$attrHash{'featID'}}; 
				$ta[8] =~ s!(ID|Name)=$s!$1=$v!ig; 
			}
		} elsif ( lc($ta[2]) eq 'mrna' ) {
			if ( defined $o2n_m{$attrHash{'featID'}} ) {
				my $s=$attrHash{'featID'}; 
				my $v=$o2n_m{$attrHash{'featID'}}; 
				$ta[8] =~ s!(ID|Name)=$s!$1=$v!ig; 
				my @aa = sort { $attrHash{'parentID'}{$a} <=> $attrHash{'parentID'}{$b} } keys %{$attrHash{'parentID'}}; 
				if (@aa == 1) {
					my $ps = $aa[0]; 
					my $pv = $o2n_g{$ps}; $pv //= $ps; 
					$ta[8] =~ s!(Parent=)$ps!$1$pv!ig; 
				} elsif ( @aa == 0 ) {
					&tsmsg("[Wrn] There is no parent for mrnaID [$s] found: $line_txt\n"); 
				} else {
					for (my $i=0; $i<$#aa; $i++) {
						my $ps = $aa[$i]; 
						my $pv = $o2n_g{$ps}; $pv //= $ps; 
						$ta[8] =~ s!(Parent=[^;]*)$ps,!$1$pv,!ig; 
					}
					my $ps = $aa[-1]; 
					my $pv = $o2n_g{$ps}; $pv //= $ps; 
					$ta[8] =~ s!(Parent=[^;]+)$ps(;|$)!$1$pv$2!ig; 
				}
			}
		} elsif ( lc($ta[2]) =~ m/^(exon|cds|five_prime_utr|three_prime_utr)$/i  ) {
			if ( $ta[8] =~ m!(?:ID|Parent)=([^\s=;]+)!i ) {
				my $s=$1; $s =~ s![\s;]+$!!; $s =~ s!:(exon|cds|five_prime_utr|three_prime_utr)(:\d+)?$!!i; 
				my $v=$o2n_m{$s}; $v //= $s; 
				if ( $s ne $v ) {
					$ta[8] =~ s!(ID|Parent|Name)=$s!$1=$v!ig; 
				}
			}
		} else {
			&stopErr("[Err] Unknown feature type [$ta[2]] in line $line_txt\n"); 
		}
		$line_txt = join("\t", @ta); 
	}
	
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'sort_by'=>$opts{'sortGffBy'} ); 
	
	return(); 
}# action_ch_makerID() 

sub action_addFaToGff {
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'seq_href'=>\%in_seq, 'sort_by'=>$opts{'sortGffBy'} ); 
}# action_addFaToGff()

sub action_compare2gffC {
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
}# action_compare2gffC() 

sub action_gffInRegion {
	# $opts{'idType'} //= 'mRNA'; 
	my $tmp_idType = (split(/,/, $opts{'idType'}))[0]; 
	my $regionFh = &openFH($opts{'gffInRegion'}, '<'); 
	my $tmpGff = &fileSunhh::new_tmp_file(); 
	my $tmpGfh = &openFH( $tmpGff, '>' ); 
	while (<$regionFh>) {
		chomp; 
		m!^\s*(#|$)! and next; 
		my @ta = split(/\t/, $_); 
		$ta[1] =~ m!^[+-]?[\d\.]+$! or do { &tsmsg("[Wrn] Skip line: $_\n"); next; }; 
		$ta[2] =~ m!^[+-]?[\d\.]+$! or do { &tsmsg("[Wrn] Skip line: $_\n"); next; }; 
		$ta[1] > $ta[2] and @ta[1,2] = @ta[2,1]; 
		my $str = $ta[3] // '.'; 
		$opts{'bothStr'} and $str = '.'; 
		# seqid / source / type / start / end / score / strand / phase / attributes(ID=???)
		if ( $str eq '.' ) {
			print {$tmpGfh} join("\t", $ta[0], 'list', $tmp_idType, $ta[1], $ta[2], '.', '+', '.', "ID=$.p")."\n"; 
			print {$tmpGfh} join("\t", $ta[0], 'list', $tmp_idType, $ta[1], $ta[2], '.', '-', '.', "ID=$.m")."\n"; 
			
		} elsif ( $str eq '+' or $str eq '-' ) {
			print {$tmpGfh} join("\t", $ta[0], 'list', $tmp_idType, $ta[1], $ta[2], '.', $str, '.', "ID=$.")."\n"; 
		} else {
			&tsmsg("[Err] Bad strand information ($str), set as both.\n"); 
			print {$tmpGfh} join("\t", $ta[0], 'list', $tmp_idType, $ta[1], $ta[2], '.', '+', '.', "ID=$.p")."\n"; 
			print {$tmpGfh} join("\t", $ta[0], 'list', $tmp_idType, $ta[1], $ta[2], '.', '-', '.', "ID=$.m")."\n"; 
		}
	}
	close($tmpGfh); 
	close ($regionFh); 
	my ( $co_gff_href ) = $gff_obj->read_gff3File('gffFile'=>$tmpGff, 'top_hier'=>$opts{'gff_top_hier'}); 
	my $rmOvlapRatio1 = -1; $opts{'onlyFull'} and $rmOvlapRatio1 = 1; 
	my $rmOvlapStrand = 'Single' ; 
	my ($kept_topIDs_aref) = &getID_rmOvlap_gff(
	 \%in_gff, $co_gff_href, 
	 'rmOvlapRatio1' => $rmOvlapRatio1, 
	 'rmOvlapType'   => $opts{'idType'}, 
	 'rmOvlapStrand' => $rmOvlapStrand, 
	 'getOvlap'      => 1 
	); 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
	
	unlink($tmpGff); 
}# action_gffInRegion() 

sub action_ovlapLongest {
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
}# action_ovlapLongest() 

sub action_islandGene {
	$opts{'islandFeatType'} //= 'CDS,exon,match_part'; 
	$opts{'islandStrand'} //= 'Both'; 
	
	my ($kept_topIDs_aref) = &getID_islandGene_gff(
	 \%in_gff, 
	 'islandFeatType' => $opts{'islandFeatType'}, 
	 'islandStrand'   => $opts{'islandStrand'}, 
	 'islandFlank' => $opts{'islandGene'}
	); 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
}# action_islandGene() 

sub action_sort {
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'sort_by'=>$opts{'sortGffBy'} ); 
}# action_sort() 
sub action_gffret {
	my $kept_topIDs_aref = []; 
	my $lisFh = &openFH($opts{'gffret'}, '<'); 
	while (<$lisFh>) {
		chomp; m!^\s*($|#)! and next; 
		my @ta = split(/\t/, $_); 
		push(@$kept_topIDs_aref, $ta[0]); 
	}
	close ($lisFh); 
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'topIDs_aref'=>$kept_topIDs_aref, 'sort_by'=>$opts{'sortGffBy'} ); 
}# action_gffret() 

sub action_listTopID {
	for ( sort { $in_gff{'ID2lineN'}{$a}<=>$in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
		my $ln = $in_gff{'ID2lineN'}{$_}; 
		my $outFeatID = $in_gff{'lineN2hash'}{ $ln }{'attrib'}{'featID'} // $_; 
		print {$oFh} "$outFeatID\n"; 
	}
}# action_listTopID () 


sub action_getLoc {
	my @chkTypes = map { lc($_) } split(/,/, $opts{'getLoc'}); 
	my %oseq_ID; 
	if ( $opts{'joinLoc'} ) {
		# featParent_ID\\tchr_ID\\tfeat_Start_min\\tfeat_End_max\\tStrand(+/-)\tfeat_Start_1,feat_End_1;feat_Start2,feat_End_2;...
		print {$oFh} join("\t", qw/LenInFeat ParentID SeqID Start End Strand Blocks/)."\n"; 
		my %ok_parentID; 
		for my $featID ( sort { $in_gff{'ID2lineN'}{$a} <=> $in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
			my $featLn = $in_gff{'ID2lineN'}{$featID}; 
			defined $in_gff{'lineN2hash'}{$featLn} or next; 
			my @parentIDs = sort { $a cmp $b } keys %{$in_gff{'CID2PID'}{$featID}}; 
			@parentIDs > 0 or @parentIDs = ( $featID ); 
			for my $parentID ( @parentIDs ) {
				$ok_parentID{$parentID} and next; 
				my %typeLocLis; 
				for my $tmp_type (@chkTypes) {
					%typeLocLis = &id_to_locLis(\%in_gff, $parentID, [$tmp_type]); 
					if (@{$typeLocLis{'locLis'}} > 0) {
						# @arr = ( $topID, $typeLocLis{'min'}, $typeLocLis{'max'}, $typeLocLis{'strand'}, $typeLocLis{'featLen'} ); 
						last; 
					}
				}
				if ( @{$typeLocLis{'locLis'}} > 0 ) {
					my @ta; 
					for my $tr ( @{$typeLocLis{'locLis'}} ) {
						push(@ta, "$tr->[0],$tr->[1]"); 
					}
					print {$oFh} join("\t", $typeLocLis{'featLen'}, $parentID, @typeLocLis{qw/scfID min max strand/}, join(";", @ta))."\n"; 

					if ( defined $oFasFh and defined $in_seq{$typeLocLis{'scfID'}} ) {
						my $tseq = ''; 
						for my $tr (@{$typeLocLis{'locLis'}}) {
							my $tseq_1 = substr( $in_seq{ $typeLocLis{'scfID'} }, $tr->[0]-1, $tr->[1]-$tr->[0]+1 ); 
							if ( $typeLocLis{'strand'} eq '+' ) {
								; 
							} elsif ( $typeLocLis{'strand'} eq '-' ) {
								&fastaSunhh::rcSeq( \$tseq_1, 'rc' ); 
							} else {
								&stopErr("[Err] Unknown strand [$typeLocLis{'strand'}] which should be '+'/'-'\n"); 
							}
							$tseq .= $tseq_1; 
						}
						my $tk = $parentID; 
						while (defined $oseq_ID{$tk}) { $tk .= "::rep"; }
						$tseq =~ s!(.{50})!$1\n!g; chomp($tseq); 
						print {$oFasFh} ">$tk\n$tseq\n"; 
					}
				} else {
					# &tsmsg("[Wrn] No information found for parentID [$parentID]\n"); 
				}
				$ok_parentID{$parentID} = 1; 
			}
		}# End for my $featID
	} else {
		# feat_ID\\tfeatParent_ID\\tSeq_ID\\tfeat_Start\\tfeat_End\\tStrand(+/-)
		print {$oFh} join("\t", qw/FeatID ParentID SeqID Start End Strand/)."\n"; 
		for my $featID ( sort { $in_gff{'ID2lineN'}{$a} <=> $in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
			my $featLn = $in_gff{'ID2lineN'}{$featID}; 
			defined $in_gff{'lineN2hash'}{$featLn} or next; 
			my $feat_href = $in_gff{'lineN2hash'}{$featLn}; 
			my $rawFeatID = $feat_href->{'attrib'}{'featID'} // $featID; 
			my @parentIDs = sort { $a cmp $b } keys %{$in_gff{'CID2PID'}{$featID}}; 
			@parentIDs > 0 or @parentIDs = ( $featID ); 
			my %typeLocLis; 
			for my $tmp_type (@chkTypes) {
				%typeLocLis = &id_to_locLis(\%in_gff, $featID, [$tmp_type]); 
				if ( @{$typeLocLis{'locLis'}} > 0 ) {
					last; 
				}
			}
			@{$typeLocLis{'locLis'}} > 0 or next; 
			for my $parentID ( @parentIDs ) {
				for my $tr ( @{$typeLocLis{'locLis'}} ) {
					print {$oFh} join( "\t", $rawFeatID, $parentID, $typeLocLis{'scfID'}, $tr->[0], $tr->[1], $tr->[2])."\n"; 
					# Edit here. 
					if ( defined $oFasFh and defined $in_seq{$typeLocLis{'scfID'}} ) {
						my $tseq = ''; 
						my $tseq_1 = substr( $in_seq{ $typeLocLis{'scfID'} }, $tr->[0]-1, $tr->[1]-$tr->[0]+1 ); 
						# Strand could be stored in $typeLocLis{'strand'} and in $tr->[2] ; 
						if ( $tr->[2] eq '+' ) {
							; 
						} elsif ( $tr->[2] eq '-' ) {
							&fastaSunhh::rcSeq( \$tseq_1, 'rc' ); 
						} else {
							&stopErr( "[Err] Unknown strand [$tr->[2]] which should be '+'/'-'\n" ); 
						} 
						$tseq = $tseq_1; 
						my $tk = $rawFeatID; 
						while ( defined $oseq_ID{$tk} ) { $tk .= "::rep"; }
						$tseq =~ s!(.{50})!$1\n!g; chomp($tseq); 
						print {$oFasFh} ">$tk\n$tseq\n"; 
					}# End if ( defined $oFasFh ... 
				}# End for my $tr ( @{$typeLocLis{'locLis'}} ) 
			}# End for my $parentID ( @parentIDs ) 
		}
	}# End if ( $opts{'joinLoc'} ) else 
}# action_getLoc() 

sub action_list_intron {
	print {$oFh} join("\t", qw/ParentID Intron_Order ChrID Start End Strand/)."\n"; 
	my %ok_parentID; 
	for my $featID ( sort { $in_gff{'ID2lineN'}{$a} <=> $in_gff{'ID2lineN'}{$b} } keys %{$in_gff{'lineN_group'}} ) {
		my $featLn = $in_gff{'ID2lineN'}{$featID}; 
		defined $in_gff{'lineN2hash'}{$featLn} or next; 
		my $feat_href = $in_gff{'lineN2hash'}{$featLn}; 
		my $rawFeatID = $feat_href->{'attrib'}{'featID'} // $featID ; 
		my @parentIDs = sort { $a cmp $b } grep { !(defined $ok_parentID{$_}) } keys %{$in_gff{'CID2PID'}{$featID}} ; 
		# my @parentIDs = sort { $a cmp $b } grep { !(defined $ok_parentID{$_}) } keys %{$in_gff{'lineN2hash'}{$featLn}{'attrib'}{'parentID'}} ; 
		@parentIDs > 0 or next; 
		# @parentIDs > 0 or @parentIDs = ( $featID ); 
		@parentIDs > 1 and do { &tsmsg("[Wrn] Be careful! There are more than 1 parents for feature [$rawFeatID]\n"); }; 
		for my $tid1 (@parentIDs) {
			my %exLocLis = &id_to_locLis(\%in_gff, $tid1, [ lc($opts{'intron_byFeat'}) ]); 
			#@{$exLocLis{'locLis'}} > 0 or next; 
			&fix_locLis( \%exLocLis ); 
			my %intronLoc = &interval_locLis( $exLocLis{'locLis'}, 0 ); 
			$exLocLis{'strand'} =~ m/^[+-]$/ or &stopErr("[Err] Bad strand [$exLocLis{'strand'}]\n"); 
			$exLocLis{'strand'} eq '-' and @{$intronLoc{'locLis'}} = reverse( @{$intronLoc{'locLis'}} ); 
			for (my $j=0; $j<@{$intronLoc{'locLis'}}; $j++) {
				my $k = $j+1; 
				print {$oFh} join("\t", 
				 $tid1, 
				 $k, 
				 $exLocLis{'scfID'}, 
				 $intronLoc{'locLis'}[$j][0], 
				 $intronLoc{'locLis'}[$j][1], 
				 $exLocLis{'strand'} 
				) . "\n"; 
			}
			$ok_parentID{$tid1} = 1; 
		}
	}
	return (); 
}# action_list_intron() 

sub action_fixTgt {
	for my $line ( grep {defined $_ and $_ !~ m!^\s*(#|$)!} values %{$in_gff{'lineN2line'}} ) {
		my @ta = split(/\t/, $line); 
		$ta[8] =~ s!(Target=\s*\S+\s+)0+(\s+\d+)!${1}1${2}!i; 
		$line = join("\t", @ta); 
	}
	$gff_obj->write_gff3File( 'outFH'=>$oFh, 'gff3_href'=>\%in_gff, 'sort_by'=>$opts{'sortGffBy'} ); 
}# action_fixTgt() 

################################################################################
############### StepX. inner sub-routines here. 
################################################################################

# Input  : (\%out_of_id_to_locLis)
# Return : () , the input hash \%out_of_id_to_locLis will be changed (sorted by 'strand'). 
sub fix_locLis {
	my ($hr) = @_; 
	for my $a1 ( @{$hr->{'locLis'}} ) {
		$a1->[0] > $a1->[1] and ($a1->[0], $a1->[1]) = ($a1->[1], $a1->[0]); 
	}
	if ( $hr->{'strand'} eq '+' ) {
		@{$hr->{'locLis'}} 
		= 
		map {
		 $_->[2] //= $hr->{'strand'}; 
		 $_->[2] eq $hr->{'strand'} or do { &tsmsg("[Wrn] Change feature_strand [$_->[2]] to [$hr->{'strand'}]\n"); $_->[2] = $hr->{'strand'}; }; 
		 $_; 
		} sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$hr->{'locLis'}}; 
	} elsif ( $hr->{'strand'} eq '-' ) {
		@{$hr->{'locLis'}} 
		= 
		map {
		 $_->[2] //= $hr->{'strand'}; 
		 $_->[2] eq $hr->{'strand'} or do { &tsmsg("[Wrn] Change feature_strand [$_->[2]] to [$hr->{'strand'}]\n"); $_->[2] = $hr->{'strand'}; }; 
		 $_; 
		} sort { $b->[1] <=> $a->[1] || $b->[0] <=> $a->[0] } @{$hr->{'locLis'}}; 
	} else {
		&stopErr("[Err] Unknown strand [$hr->{'strand'}] [$hr->{'locLis'}[0] $hr->{'locLis'}[1]]\n"); 
	}
	return(); 
}# fix_locLis() 


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


# Input : (\%in_gff_1, \%in_gff_2_Compare, 'rmOvlapLen'=>-1, 'rmOvlapRatio1'=>-1, 'rmOvlapRatio2'=>-1, 'rmOvlapType'=>'exon,CDS,match_part', 'rmOvlapStrand'=>'Both', 'getOvlap'=>0)
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
	$parm{'getOvlap'} //= 0; 
	
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
	
	if ( $parm{'getOvlap'} ) {
		my %rmID = map { $_ => 1 } @back_topIDs; 
		my @useIDs = grep { !(defined $rmID{$_}) } @topIDs_1; 
		@back_topIDs = @useIDs; 
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
	# push( @off_lines_hash, map { $gff->{'lineN2hash'}{$_} } grep { defined $gff->{'lineN2hash'}{$_} } @{ $gff->{'lineN_group'}{$featID}{'curLn'} } ); 
	# push( @off_lines_hash, map { $gff->{'lineN2hash'}{$_} } grep { defined $gff->{'lineN2hash'}{$_} } @{ $gff->{'lineN_group'}{$featID}{'offLn'} } ); 
	defined $gff->{'lineN2hash'}{$curLn} and push( @off_lines_hash, $gff->{'lineN2hash'}{$curLn} ); 
	push( @off_lines_hash, 
	 map { $gff->{'lineN2hash'}{$_} } grep { defined $gff->{'lineN2hash'}{$_} } sort { $a<=>$b } map { $gff->{'ID2lineN'}{$_} } keys %{ $gff->{'PID2CID'}{$featID} }
	); 
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

sub load_gff_fas {
	my ($gff_href, $seq_href) = @_; 
	$gff_href //= {}; 
	$seq_href //= {}; 
	my %in_gff = %$gff_href; 
	my %in_seq = %$seq_href; 
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
	%$gff_href = %in_gff; 
	%$seq_href = %in_seq; 
	return ($gff_href, $seq_href); 
}# load_gff_fas

######################################################################
######################################################################
######################################################################
######################################################################
#!/usr/bin/perl


sub rev_comp {
	my @back; 
	for (@_) {
		my $ts = reverse($_); 
		$ts =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
		push(@back, $ts); 
	}
	return @back; 
}
