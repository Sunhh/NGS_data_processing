package mcsSunhh; 

use strict; 
use warnings; 
our @EXPORT = qw(); 
our @EXPORT_OK = qw(); 

use mathSunhh; 
use fileSunhh; 
use LogInforSunhh; 

############################################################
#  Methods
############################################################

=head1 chrP_to_plotP( 'chrLis' => &mcsSunhh::_readInChrLis(), 'chrID' => $chr_ID , 'chrP' => $chr_posi <1> , 'beginPlotP' => $zeroPlotV <0> , 'BpPoint' => $bp_per_point , 'SepBorder' => $dist <0> )

Input      : Required - chrLis , chrID 
  beginPlotP is the plot-position for zero chr-position of the first chrID . 
  SepBorder is the plot-length between two plotted chromosomes. 

Return     : ( $plotPosi_1, $plotPosi_2, ... )

  If there is no $chrID found , return an empty array. 

=cut
sub chrP_to_plotP {
	my %parm = &mathSunhh::_setHashFromArr( @_ ); 
	$parm{'SepBorder'} //= 0; 
	$parm{'BpPoint'} //= 1; 
	defined $parm{'chrLis'} or &stopErr("[Err] I need 'chrLis' from mcsSunhh::_readInChrLis()\n"); 
	my @backP; 
	( defined $parm{'chrID'} and defined $parm{'chrLis'}{'ord'}{$parm{'chrID'}} ) or return(@backP); 
	$parm{'chrP'} //= 1; 
	$parm{'beginPlotP'} //= 0; 
	my $prevP = $parm{'beginPlotP'}; 
	for my $t_chr (@{ $parm{'chrLis'}{'arr'} }) {
		if ( $t_chr->[0] eq $parm{'chrID'} ) {
			push( @backP, $prevP + $parm{'chrP'}/$parm{'BpPoint'} ); 
		}
		$prevP += ($t_chr->[1] / $parm{'BpPoint'}); 
		$prevP += $parm{'SepBorder'}; 
	}
	
	return(@backP); 
}# chrP_to_plotP


=head1 _readInMcsGff( $fileName_mcs_gff ) 

Format of $fileName_mcs_gff : ChrID \\t GeneID \\t Start \\t End
This file should be already sorted. 

Return    : (\%chr2gen, \%gen2loc)
 %chr2gen : {chromID} => [ [genID, chrS, chrE], [], ... ]
 %gen2loc : {genID}   => [chrID, chrS, chrE]
=cut
sub _readInMcsGff {
	my $inFh = &openFH(shift, '<'); 
	my (%chr2gen, %gen2loc); 
	while ( &wantLineC($inFh) ) {
		my ( $chrID, $genID, $chrS, $chrE ) = &splitL("\t", $_); 
		push( @{$chr2gen{$chrID}}, [$genID, $chrS, $chrE] ); 
		defined $gen2loc{$genID} and &tsmsg("[Msg] genID [$genID] repeats, and loc infor masked by latter one\n"); 
		$gen2loc{$genID} = [ $chrID, $chrS, $chrE ]; 
	}
	close($inFh); 
	return (\%chr2gen, \%gen2loc); 
}# _readInMcsGff() 


=head1 _update_alnInfoText_byPair( $alnInfo[$idx_num]{'pair'}, $idx_num, $alnInfo{$idx_num}{'info'} )

Input    : A element of @alnInfo from &_readInAln(), which is not the header part. 

  $alnInfo{$idx_num}{'info'} is not necessary, but should not be given wrongly. 
  If providing {'info'} [alnID, score, evalue, Num_pairs, chrID_1, chrID_2, strand(plus|minus)], 
    the values could be [-1,   999999, 0,      -1,        XX,      XX,      plus], 
    '-1' means undefined. 

Return   : ( $text )
  A new generated $alnInfo{$idx_num}{'text'}

=cut
sub _generate_alnInfoText_byPair {
	my ( $pAR, $alnID, $iAR ) = @_; 
	my $n_pairs = scalar(@$pAR); 
	$iAR->[0] //= $alnID; 
	$iAR->[0] == -1 and $iAR->[0] = $alnID; 
	$iAR->[0] == $alnID or &stopErr("[Err] alnID [$alnID] is different from given iAR [$iAR->[0]]\n"); 
	$iAR->[1] //= '999999'; # score
	$iAR->[2] //= '0'; # e_value
	$iAR->[3] //= $n_pairs; # N
	$iAR->[3] == -1 and $iAR->[3] = $n_pairs; 
	$iAR->[3] == $n_pairs or &stopErr("[Err] Num_pairs [$n_pairs] is different from given iAR [$iAR->[3]]\n"); 
	$iAR->[4] //= 'XX'; # chr1
	$iAR->[5] //= 'XX'; # chr2
	$iAR->[6] //= 'plus'; # plus|minus
	my $text = ''; 
	$text .= "## Alignment $iAR->[0]: score=$iAR->[1] e_value=$iAR->[2] N=$iAR->[3] $iAR->[4]\&$iAR->[5] $iAR->[6]\n"; # ## Alignment 0: score=7185.0 e_value=0 N=148 Cma_Chr04&Cma_Chr04 plus
	# Check if we have ka/ks information in gene pairs; 
	my $has_ks = 0; 
	for (@$pAR) {
		defined $_->[3] and $_->[3] ne '' and do { $has_ks = 1; last; }; 
		defined $_->[4] and $_->[4] ne '' and do { $has_ks = 1; last; }; 
		defined $_->[5] and $_->[5] ne '' and do { $has_ks = 1; last; }; 
	}
	for (my $i=0; $i<$n_pairs; $i++) {
		my $j=$i+1; 
		my @tp = @{$pAR->[$i]}; 
		if ($has_ks) {
			$text .= join( "\t", sprintf("%3d-%3d:", $iAR->[0], $j), @tp[0 .. 2], $tp[4], $tp[3] ) . "\n"; # geneID_1, geneID_2, pair_evalue, Ks, Ka ; 
		} else {
			$text .= join( "\t", sprintf("%3d-%3d:", $iAR->[0], $j), @tp[0 .. 2] ) . "\n"; # geneID_1, geneID_2, pair_evalue
		}
	}
	return ($text); 
} # _update_alnInfoText_byPair () 

=head1 _readInAln( $fileName_mcs_aln, $useYN <1> )

Input    : Output .collinearity file of MCScanX ; 
           If $useYN is FALSE (0) and Nei-Gojobori (NG) Ks values exist, I will try to use NG_Ks instead of default Yang-Nielson (YN) Ks values. 

Return   : ( \@alnInfo )
  @alnInfo : 
    [idx_num]{'info'} => [alnID, score, evalue, Num_pairs, chrID_1, chrID_2, strand(plus|minus)]
    [idx_num]{'pair'} => [ [genID_1, genID_2, pair_evalue, Ka, Ks, Ka/Ks], [], ... ]
    [idx_num]{'text'} => $text_of_current_block

=cut
sub _readInAln{
	my ($fn, $useYN) = @_; 
	$useYN //= 1; 
	my $inFh = &openFH($fn, '<'); 
	my @alnInfo; 
	my $aln_idx = 0; 
	while (<$inFh>) {
		chomp; 
		if ( m!^\s*$|^####|^# \S+:\s+\S+$|# (MAX GAPS|Number of collinear genes|Number of all genes):! ) {
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
			next; 
		}
		if ( m!^## Alignment ! ) {
			# ## Alignment 0: score=7137.0 e_value=0 N=147 Ma1&Ma1 plus
			m!## Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) ([^\&\s]+)\&([^\&\s]+) (plus|minus|X+)! or die "ALN:$_\n"; 
			my ($alnID, $score, $eval, $n, $chr1, $chr2, $str)
			=  ($1,     $2,     $3,    $4, $5,    $6,    $7); 
			$aln_idx++; 
			$alnInfo[$aln_idx]{'info'} = [$alnID, $score, $eval, $n, $chr1, $chr2, $str]; 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} elsif ( m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+))?)?$! ) { 
			# #  0-  0:        Cma_000007      Cma_000973        2e-57 yn_Ka yn_Ks yn_kaks
			$aln_idx > 0 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tka, $tks, $tw) 
			= 
			   ($1,     $2,        $3,    $4,    $5,    $6,   $7,   $8); 
			$tka //= ''; $tks //= ''; $tw //= ''; 
			$alnID == $alnInfo[$aln_idx]{'info'}[0] or &stopErr("[Err] line_alnID=$alnID not fitting upper level (alnID=$alnInfo[$aln_idx]{'info'}[0]).\n"); 
			push(@{$alnInfo[$aln_idx]{'pair'}}, [$gid1, $gid2, $eval, $tka, $tks, $tw]); 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} elsif ( m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+)\t(\S+))?)?$! ) { 
			# #  0-  0:  Cma_000007      Cma_000973        2e-57 0.5258  0.0154  0.5387  0.0148
			#                                                    yn_Ks   yn_Ka   ng_Ks   ng_Ka
			# This is to fit result from haibao tang's python ks calculation. 
			$aln_idx > 0 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tks, $tka, $t_ngKs, $t_ngKa) 
			= 
			(   $1,     $2,        $3,    $4,    $5,    $6,   $7,   $8,      $9); 
			my $tw; 
			if ( !$useYN and defined $t_ngKs ) {
				$tks = $t_ngKs; 
				$tka = $t_ngKa; 
			}
			$tks < 0 and $tks = 'nan'; 
			if ( defined $tks and $tks ne 'nan') {
				$tw = ( $tks > 0 ) ? $tka/$tks : 'nan'; 
			}
			$tks eq 'nan' and $tw = 'nan'; 
			$tka //= ''; $tks //= ''; $tw //= ''; 
			$alnID == $alnInfo[$aln_idx]{'info'}[0] or &stopErr("[Err] line_alnID=$alnID not fitting upper level (alnID=$alnInfo[$aln_idx]{'info'}[0]).\n"); 
			push(@{$alnInfo[$aln_idx]{'pair'}}, [$gid1, $gid2, $eval, $tka, $tks, $tw]); 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} else {
			&stopErr("[Err] Unable to parse line: $_\n"); 
		}
	}
	close($inFh); 
	return(\@alnInfo); 
}# _readInAln


=head1 _readInAlnTbl( $fileName_aln2table )

Return       : ( \@tabInfo )

  @tabInfo = ( \@blk_1, \@blk_2, ... )
    @blk_1 = ( $blkID , 
      [$chrID_1, $chrS_1, $chrE_1], 
      [$chrID_2, $chrS_2, $chrE_2], 
      $strand, $alnScore, $alnEvalue, $alnNumber, 
      [$gene1_1, $gene1_2, ...], 
      [$gene2_1, $gene2_2, ...], 
      [$Ka_1, $Ka_2, ...], 
      [$Ks_1, $Ks_2, ...], 
      [$KaKs_1, $KaKs_2, ...], 
      all of the rest
    )
   If $blkID == 'BlkID', this block should be the header of the file. 
=cut
sub _readInAlnTbl {
	my ($fn) = @_; 
	my $fh = &openFH($fn, '<'); 
	my @back; 
	while ( &wantLineC($fh) ) {
		my @ta = &splitL("\t", $_); 
		my @tb; 
		$tb[0] = $ta[0]; 
		$tb[1] = [@ta[1,2,3]]; 
		$tb[2] = [@ta[4,5,6]]; 
		@tb[3,4,5,6] = @ta[7,8,9,10]; 
		$tb[7] = [ &splitL(",", $ta[11]) ]; 
		$tb[8] = [ &splitL(",", $ta[12]) ]; 
		$tb[9] = [ &splitL(",", $ta[13]) ]; 
		$tb[10] = [ &splitL(",", $ta[14]) ]; 
		$tb[11] = [ &splitL(",", $ta[15]) ]; 
		push(@tb, @ta[16 .. $#ta]); 
		push(@back, [@tb]); 
	}
	close($fh); 
	return(\@back); 
}# _readInAlnTbl()

=head1 _readInChrLis ( $file_name )

Input        : Input file format : Chr_ID \\t Chr_Len \\n

Return       : ( \%back )
  'arr' => ([ [chr1_ID, chr1_Len, cum_len_prevE, repeat_Idx], [chr2_ID, chr2_Len, cum_len_prevE, repeat_Idx], ... ]); 
  'ord' => { chr_ID => [idx_in_arr, ...] }
	'repN'=> { chr_ID => $repeat_time }

=cut
sub _readInChrLis {
	# $_[0] : input file name; Format: Scf_ID \\t  Scf_Len  
	my $inFh = &openFH($_[0], '<'); 
	my %back; 
	my $cum_len = 0; 
	while (<$inFh>) {
		chomp; 
		m/^\s*(#|$)/ and next; 
		my @ta = split(/\t/, $_); 
		( defined $ta[1] and $ta[1] > 0 ) or do { &tsmsg("[Wrn] Skip scaffold [$ta[0]] with bad length [$ta[1]].\n"); next; }; 
		$back{'repN'}{$ta[0]}++; 
		defined $back{'ord'}{$ta[0]} and &tsmsg("[Wrn] Repeated scaffold ID [$ta[0]]\n"); 
		push(@{$back{'arr'}}, [ $ta[0], $ta[1], $cum_len, $back{'repN'}{$ta[0]}-1 ]); 
		push(@{ $back{'ord'}{$ta[0]} }, $#{ $back{'arr'} } ); 
		$cum_len += $ta[1]; 
	}
	close($inFh); 
	return(\%back); 
}# _readInChrLis() 


1; 