#!/usr/bin/perl
# 2016-03-23 This is used to count read number within genes for rnaseq programs. 
#   The motivation to re-write this function is to shorten processing time. 
#   Because time is limit, and because my rnaseq projects are usually strand-specific and single stranded, here I only write counting for that type. 
#   I may add function for other strand-specific/nonspecific counting in the future. 
use strict; 
use warnings; 
use SeqAlnSunhh; 
use mathSunhh; 
use fileSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 

	"outPref:s",        # 'out'. Out prefix. 

	"highMemory!",      # To trigger this on will cause higher memory usage. 
	"onlyCntRdSE!",     # Only count a read if this read's start/end position locates in a gene. 
	"withYiBinMethod!", # Use Yi Zheng's bin method which is only suitable for -onlyCntRdSE counting. 
	"senseStrand:s",    # Default 'R'; 'R' - for different strand of gene , 'F' - for same strand of gene. 
	"OnlyCntTotal!",    # Only count total reads discarding if it locates in a gene. 

	# converters 
	"loc_1to1:s",       # filename. file format : old_ID \\t old_posi \\t new_ID \\t new_posi
	"bam_isList!", 

	# Filters 
	"max_mismatchN:i",  # -1 . could be [0-...]
	"max_mismatchR:f",  # -1 . could be [0-1]

	"exe_samtools:s",     # samtools_1.3 
); 

$opts{'withYiBinMethod'} and $opts{'onlyCntRdSE'} = 1; 
$opts{'senseStrand'}   //= 'R'; 
$opts{'exe_samtools'}  //= 'samtools'; 
$opts{'max_mismatchN'} //= -1; 
$opts{'max_mismatchR'} //= -1; 

$opts{'outPref'} //= 'out'; 

my $help_txt = <<HH; 

perl $0 in.bed_fit_YZpipe in.bam   1>in.bam.cnt   2>in.bam.cnt.err

-outPref                    ['out'] Will output 'outPref'.sense.cnt and 'outPref'.anti.cnt files. 

-highMemory                 [Boolean]

-withYiBinMethod            [Boolean]
-onlyCntRdSE                [Boolean]

-senseStrand                ['R'] R/F

-loc_1to1                   [filename] file format : old_ID \\t old_posi \\t new_ID \\t new_posi
-bam_isList                 [Boolean] in.bam is a bam list [ bam_filename \\t out_prefix ] if given. 

-max_mismatchN              [-1]
-max_mismatchR              [-1]

-exe_samtools               ['samtools']

-help                       [Boolean]

HH

@ARGV < 2 and &LogInforSunhh::usage( $help_txt ); 
$opts{'help'} and &LogInforSunhh::usage( $help_txt ); 

# [Sunhh@Falcon temp]$ head -3 r6_maker_final_noFa.fit_KTools.bed
# S800001_pilon   4845    5136    Cma_000001      291     -
# S800001_pilon   7100    13333   Cma_000002      3579    +
# S800001_pilon   14478   19291   Cma_000003      879     +
# [Sunhh@Falcon temp]$ samtools_1.3 view flowers_rep1.fq_hisat.srt.bam | head -3
# HWI-D00656:36:C6AFMANXX:4:1206:12397:9131       0       P1Cma_Chr01     1425    255     101M    *       0       0       ATAAGACCTCACAACCTCAAATTTTTCAACGAAATTGTAAGAAACGTTCAAATTTGATGAATTCGACACATAAGTTGTTGAAGAAAAGAAGGAAAATAAAA  BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF   AS:i:0XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:101        YT:Z:UU XS:A:-  NH:i:1
# HWI-D00656:36:C6AFMANXX:4:1312:3096:4163        0       P1Cma_Chr01     1559    255     101M    *       0       0       GTTCACTGTTTTCAGCCTCGAAAATAAAAGTTGAGGCAGGCAAGATTCAAGTTGAAGTAGGCGATTGCAGTGGGGATTGAAGAAAATAAAGGTGATGAAAC  BBBBBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF   AS:i:0XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:101        YT:Z:UU XS:A:-  NH:i:1
# HWI-D00656:36:C6AFMANXX:4:1308:20875:59522      0       P1Cma_Chr01     1597    255     101M    *       0       0       GGCAAGATTCAAGTTGAAGTAGGCGATTGCAGTGGGGATTGAAGAAAATAAAGGTGATGAAACACAAATTTTAGTGTGGGGAGAAGTCAACTCCATCATCA  BBBBBFFFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFB   AS:i:0XN:i:0   XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:101        YT:Z:UU XS:A:-  NH:i:1

my $fn_bed = shift; 
my $fn_bam = shift; 

my %loc_1to1; 
defined $opts{'loc_1to1'} and &load_loc_1to1( $opts{'loc_1to1'}, \%loc_1to1 ); 

my %bed_info = %{ &load_bed($fn_bed) }; 
# Return : ( { "$chrID"=>[ [ chrS_1idx, chrE_1idx, geneID, geneLen, strand, inputOrder, $chrID ], [], ...], ""=>, ...   } )
&tsmsg("[Msg] Indexing.\n"); 
my %bed_db = %{ &index_bed( \%bed_info ) }; 

my @bam_files; 
if ($opts{'bam_isList'}) {
	my $th = &openFH($fn_bam, '<'); 
	while (&wantLineC($th)) {
		my @ta = &splitL("\t", $_); 
		( defined $ta[0] and $ta[0] ne '' ) or next; 
		( defined $ta[1] and $ta[1] ne '' ) or &stopErr("[Err] bam_list needs two colums [bam_fn, out_prefix]\n"); 
		push(@bam_files, [$ta[0], $ta[1]]); 
	}
	close($th); 
} else {
	push(@bam_files, [$fn_bam, $opts{'outPref'}]); 
}

my %flag_aln_F = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=0,4=0' ) }; 
my %flag_aln_R = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=0,4=1' ) }; 
my %flag_aln   = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=0' ) }; 
my %flag_to_str; 
for ( keys %{ &SeqAlnSunhh::mk_flag( 'keep' => '4=0' ) } ) {
	$flag_to_str{$_} = '+'; 
}
for ( keys %{ &SeqAlnSunhh::mk_flag( 'keep' => '4=1' ) } ) {
	$flag_to_str{$_} = '-'; 
}

for my $cur (@bam_files) {
	my ($cur_bam, $cur_opref) = @$cur; 

my $sam_fh = &SeqAlnSunhh::openSam( $cur_bam, undef(), { 'wiH'=>0, 'verbose'=>1, 'exe_samtools'=>$opts{'exe_samtools'} } ); 

my %cnt; 
$cnt{'log_section'} = { 'cntN_base' => 0, 'cntN_step' => 1e5 }; 
SAM_LINE:
while (<$sam_fh>) {
	&fileSunhh::log_section( $. , $cnt{'log_section'} ) and &tsmsg("[Msg] Reading $. line.\n"); 
	chomp; 
	my @ta = &splitL("\t", $_); 
	my ($aln_rdID, $aln_flag, $aln_chrID, $aln_chrP, $aln_cigar) = @ta[0,1,2,3,5]; 
	defined $flag_aln{ $aln_flag } or next SAM_LINE; 
	if ( !defined $opts{'loc_1to1'} ) {
		$opts{'OnlyCntTotal'} or defined $bed_db{$aln_chrID} or next SAM_LINE; 
		$opts{'OnlyCntTotal'} and defined $cnt{'rdID_hash'}{$aln_rdID} and next SAM_LINE; 
	}
	$cnt{'tmp_mismatchN'} = -2; 
	if ( $opts{'max_mismatchN'} >= 0 ) {
		($cnt{'tmp_mismatchN'}, $cnt{'tmp_cigarH'}) = &SeqAlnSunhh::cnt_sam_mismatch( \@ta, 'set_rna' ); 
		$cnt{'tmp_mismatchN'} <= $opts{'max_mismatchN'} or next SAM_LINE; 
	}
	if ( $opts{'max_mismatchR'} >= 0 ) {
		$cnt{'tmp_mismatchN'} == -2 and ($cnt{'tmp_mismatchN'}, $cnt{'tmp_cigarH'}) = &SeqAlnSunhh::cnt_sam_mismatch( \@ta, 'set_rna' ); 
		$cnt{'tmp_mismatchN'} <= $opts{'max_mismatchR'} * $cnt{'tmp_cigarH'}{'RdLen'} or next SAM_LINE; 
	}
	$cnt{'rdID_hash'}{$aln_rdID} ++; 
	$opts{'OnlyCntTotal'} and next SAM_LINE; 
	my %cigar_h = %{ &SeqAlnSunhh::parseCigar( $aln_cigar ) }; 
	my ($spanS, $spanE) = ( $aln_chrP, $aln_chrP+$cigar_h{'SpanRefLen'}-1 ); 
	my $spanStr = $flag_to_str{ $aln_flag }; 

	if ( defined $opts{'loc_1to1'} ) {
		my ( $s_chrID, $s_chrP, $s_str ) = &get_newLoc( \%loc_1to1, $aln_chrID, $spanS, $spanStr ); 
		my ( $e_chrID, $e_chrP, $e_str ) = &get_newLoc( \%loc_1to1, $aln_chrID, $spanE, $spanStr ); 
		if ( $s_chrID eq $e_chrID and $s_chrID ne '' ) {
			# If $s_str ne $e_str, I will use $s_str only. 
			($aln_chrID, $spanS, $spanE, $spanStr) = ( $s_chrID, $s_chrP, $e_chrP, $s_str ); 
			$spanS > $spanE and ($spanS, $spanE) = ($spanE, $spanS); 
		} else {
			next SAM_LINE; 
		}
	}

	if ($opts{'highMemory'}) {
		$cnt{'t_pos'} = "$aln_chrID:$spanS-$spanE"; 
		if ( defined $cnt{'pos2genes'}{ $cnt{'t_pos'} } ) {
			if ( $spanStr eq '+' ) {
				for my $ar1 ( @{$cnt{'pos2genes'}{ $cnt{'t_pos'} }} ) {
					defined $cnt{'geneCnt'}{ $ar1->[2] }[2]{ $aln_rdID } or $cnt{'geneCnt'}{ $ar1->[2] }[0] ++; # Fwd
					$cnt{'geneCnt'}{ $ar1->[2] }[2]{ $aln_rdID } ++; 
				}
			} elsif ( $spanStr eq '-' ) {
				for my $ar1 ( @{$cnt{'pos2genes'}{ $cnt{'t_pos'} }} ) {
					defined $cnt{'geneCnt'}{ $ar1->[2] }[3]{ $aln_rdID } or $cnt{'geneCnt'}{ $ar1->[2] }[1] ++; # Rev 
					$cnt{'geneCnt'}{ $ar1->[2] }[3]{ $aln_rdID } ++; 
				}
			}
			next SAM_LINE; 
		}
	}
	my @loc_realIdx; 
	if ( $opts{'withYiBinMethod'} ) {
		my %used; 
		@loc_realIdx = grep { $used{$_}=1; } @{ &idxPos_byYi( $bed_db{$aln_chrID}, $spanS ) }; 
		push(@loc_realIdx, grep { !(defined $used{$_}) } @{ &idxPos_byYi( $bed_db{$aln_chrID}, $spanE ) } ); 
	} elsif ( $opts{'onlyCntRdSE'} ) {
		my %used; 
		@loc_realIdx = grep { $used{$_}=1; } @{ &mathSunhh::_hasPos_inLocDb( $bed_db{$aln_chrID}, $spanS ) }; 
		push(@loc_realIdx, grep { !(defined $used{$_}) } @{ &mathSunhh::_hasPos_inLocDb( $bed_db{$aln_chrID}, $spanE ) } ); 
	} else {
		@loc_realIdx = @{ &mathSunhh::_map_loc_to_realIdx( $bed_db{$aln_chrID}, [$spanS, $spanE] ) }; 
	}
	$opts{'highMemory'} and $cnt{'pos2genes'}{ $cnt{'t_pos'} } = [ @{ $bed_info{$aln_chrID} }[@loc_realIdx] ]; 
	if ( $spanStr eq '+' ) {
		for my $ar1 (@{ $bed_info{$aln_chrID} }[@loc_realIdx]) {
			defined $cnt{'geneCnt'}{ $ar1->[2] }[2]{ $aln_rdID } or $cnt{'geneCnt'}{ $ar1->[2] }[0] ++; # Fwd
			$cnt{'geneCnt'}{ $ar1->[2] }[2]{ $aln_rdID } ++; 
		}
	} elsif ( $spanStr eq '-' ) {
		for my $ar1 (@{ $bed_info{$aln_chrID} }[@loc_realIdx]) {
			defined $cnt{'geneCnt'}{ $ar1->[2] }[3]{ $aln_rdID } or $cnt{'geneCnt'}{ $ar1->[2] }[1] ++; # Rev
			$cnt{'geneCnt'}{ $ar1->[2] }[3]{ $aln_rdID } ++; 
		}
	} else {
		next SAM_LINE; 
	}
}
close($sam_fh); 

$cnt{'sum_cnt_all'} = scalar( keys %{$cnt{'rdID_hash'}} ); 
&tsmsg("[Rec] Total $cnt{'sum_cnt_all'} reads accepted in $cur_bam\n"); 

if ( $opts{'senseStrand'} eq 'R' ) {
	for my $ar1 ( map { @{ $bed_info{$_} } } keys %bed_info ) {
		my ($chrS, $chrE, $geneID, $geneLen, $strand, $inOrder, $chrID) = @$ar1; 
		$cnt{'geneCnt'}{ $geneID }[0] //= 0; # Fwd
		$cnt{'geneCnt'}{ $geneID }[1] //= 0; # Rev
		$cnt{'geneCnt'}{ $geneID }[2] //= {}; 
		$cnt{'geneCnt'}{ $geneID }[3] //= {}; 
		$strand eq '+' and do { @{ $cnt{'geneCnt'}{ $geneID } }[ 0,1, 2,3 ] = @{ $cnt{'geneCnt'}{ $geneID } }[ 1,0, 3,2 ]; }; 
		for (keys %{ $cnt{'geneCnt'}{ $geneID }[2] }) { $cnt{'sum_cnt_sense'}{$_} = 1; }
		for (keys %{ $cnt{'geneCnt'}{ $geneID }[3] }) { $cnt{'sum_cnt_anti'}{$_}  = 1; }
	}
} elsif ( $opts{'senseStrand'} eq 'F' ) {
	for my $ar1 ( map { @{ $bed_info{$_} } } keys %bed_info ) {
		my ($chrS, $chrE, $geneID, $geneLen, $strand, $inOrder, $chrID) = @$ar1; 
		$cnt{'geneCnt'}{ $geneID }[0] //= 0; # Fwd
		$cnt{'geneCnt'}{ $geneID }[1] //= 0; # Rev
		$cnt{'geneCnt'}{ $geneID }[2] //= {}; 
		$cnt{'geneCnt'}{ $geneID }[3] //= {}; 
		$strand eq '-' and @{ $cnt{'geneCnt'}{ $geneID } }[ 0,1, 2,3 ] = @{ $cnt{'geneCnt'}{ $geneID } }[ 1,0, 3,2 ]; 
		for (keys %{ $cnt{'geneCnt'}{ $geneID }[2] }) { $cnt{'sum_cnt_sense'}{$_} = 1; }
		for (keys %{ $cnt{'geneCnt'}{ $geneID }[3] }) { $cnt{'sum_cnt_anti'}{$_}  = 1; }
	}
} else {
	&stopErr("[Err] Unknown -senseStrand type [$opts{'senseStrand'}].\n"); 
}

open O1,'>',"$cur_opref.sense.cnt" or &stopErr("[Err] Failed to write file [$cur_opref.sense.cnt]\n"); 
open O2,'>',"$cur_opref.anti.cnt" or &stopErr("[Err] Failed to write file [$cur_opref.sense.cnt]\n"); 

print O1 join("\t", 'Gene_ID',  "CntSens_$cur_bam")."\n"; 
print O1 join("\t", 'SumSense', scalar(keys %{$cnt{'sum_cnt_sense'}}))."\n"; 
print O2 join("\t", 'Gene_ID',  "CntAnti_$cur_bam")."\n"; 
print O2 join("\t", 'SumAnti', scalar(keys %{$cnt{'sum_cnt_anti'}}))."\n"; 

for my $ar1 (sort {$a->[5]<=>$b->[5]} map { @{ $bed_info{$_} } } sort keys %bed_info ) {
	my ($chrS, $chrE, $geneID, $geneLen, $strand, $inOrder, $chrID) = @$ar1; 
	print O1 join("\t", $geneID, $cnt{'geneCnt'}{ $geneID }[0])."\n"; 
	print O2 join("\t", $geneID, $cnt{'geneCnt'}{ $geneID }[1])."\n"; 
}

close(O1); 
close(O2); 

}# End for (@bam_files)

&tsmsg("[Rec] $0 done.\n"); 


sub index_bed {
	my ($bed_hr) = @_; 
	my %back; 
	for my $k1 (keys %$bed_hr) {
		&tsmsg("[Msg] Indexing $k1\n"); 
		if ( $opts{'withYiBinMethod'} ) {
			my @seLoc = map { [ $_->[0], $_->[1] ] } @{$bed_hr->{$k1}}; 
			$back{$k1}{'seLoc'} = \@seLoc; 
			$back{$k1}{'binSize'} = 1000; 
			for (my $j=0; $j<@{$back{$k1}{'seLoc'}}; $j++) {
				my $ar1 = $back{$k1}{'seLoc'}[$j]; 
				for (my $i=int($ar1->[0]/$back{$k1}{'binSize'}); $i<=int($ar1->[1]/$back{$k1}{'binSize'}); $i++) {
					push(@{$back{$k1}{'binIdx'}{$i}}, $j); 
				}
			}
		} else {
			my @seLoc = map { [ $_->[0], $_->[1] ] } @{$bed_hr->{$k1}}; 
			$back{$k1} = &mathSunhh::index_SEloc( \@seLoc ); 
		}
	}

	return(\%back); 
}# index_bed ()

sub idxPos_byYi {
	my ($ah, $posi) = @_; 
	my $binIdx = int( $posi/$ah->{'binSize'} ); 
	defined $ah->{'binIdx'}{$binIdx} or return([]); 
	my @back; 
	for my $realIdx (@{$ah->{'binIdx'}{$binIdx}}) {
		$posi >= $ah->{'seLoc'}[$realIdx][0] and $posi <= $ah->{'seLoc'}[$realIdx][1] and push(@back, $realIdx); 
	}
	return(\@back); 
}# idxPos_byYi () 


# Return : ( { "$chrID"=>[ [ chrS_1idx, chrE_1idx, geneID, geneLen, strand, inputOrder, $chrID ], [], ...], ""=>, ...   } )
sub load_bed {
	my $fh = &openFH($_[0], '<'); 
	my $cnt = -1; 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		$ta[1] ++; 
		$cnt ++; 
		push(@{$bed_info{$ta[0]}}, [ @ta[1 .. 5], $cnt, $ta[0] ]); 
	}
	close($fh); 
	return(\%bed_info); 
}# load_bed 
sub load_loc_1to1 {
	my ($fn, $hr) = @_; 
	$hr //= {}; 
	my $fh = &openFH( $fn, '<' ); 
	&tsmsg("[Msg] Loading loc_1to1 [$fn]\n"); 
	my %c; 
	$c{'log_section'} = { 'cntN_base' => 0, 'cntN_step' => 1e5 }; 
	while (&fileSunhh::wantLineC( $fh )) {
		&fileSunhh::log_section( $. , $c{'log_section'} ) and &tsmsg("[Msg] Reading $. line.\n");
		my ( $oID, $oP, $nID, $nP, $nStr ) = split(/\t/, $_); 
		push( @{$hr->{$oID}{$oP}}, [ $nID, $nP, $nStr ] ); 
	}
	close($fh); 
	&tsmsg("[Msg] Loaded [$fn]\n"); 
	return(\$hr); 
}# load_loc_1to1() 

sub get_newLoc {
	my ($hr, $oID, $oP, $oStr) = @_; # input Old coordinates. 
	$oStr =~ m/^[+-]$/ or &stopErr("[Err] Bad strand input of get_newLoc( $hr, $oID, $oP, $oStr ).\n"); 
	my ($nID, $nP, $nStr) = ('', '', ''); 
	my $flank_len = 50; 
	defined $hr->{$oID} or return($nID, $nP); 
	if (defined $hr->{$oID}{$nP}) {
		($nID, $nP, $nStr) = @{ $hr->{$oID}{$nP}[0] }; 
	} else {
		for (my $i=1; $i<=$flank_len; $i++) {
			my $tP_1 = $oP-$i; 
			defined $hr->{$oID}{$tP_1} and do { ($nID, $nP, $nStr) = @{ $hr->{$oID}{$tP_1}[0] }; last; }; 
			my $tP_2 = $oP+$i; 
			defined $hr->{$oID}{$tP_2} and do { ($nID, $nP, $nStr) = @{ $hr->{$oID}{$tP_2}[0] }; last; }; 
		}
	}
	$oStr eq '-' and $nStr =~ tr/+-/-+/; 
	return($nID, $nP, $nStr); 
}# get_newLoc() 

