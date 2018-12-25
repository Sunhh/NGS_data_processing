#!/usr/bin/perl
# 2018-10-26 : step 3 . Check if possible transmitted reads are contained by non/self-grafted samples. 
# 2018-10-30 : -transRdMaxLen cut transmitted read to shorter length to fit -bgRdBam ; 
#   Input : -transRdBam:s , -bgRdBam:s@ , -pref , -tgt_fa , -exe_samtools , -wrk_dir , -outBam ; 
#            pref.src2tgt_rd.bam ; pref_comb.nonTransRmdup.bam ; 
#   Output: wrk_dir/pref.src2tgt_cleanRd.bam
# Step 1 : Remove RG/PG information from transRdBam, and add RG:Z:src; 
#          Remove RG/PG information from bgRdBam, and add RG:Z:tgt; 
#          Combine transRdBam and bgRdBam, and then sort them by position; 
# Step 2 : Check RG_src_reads one by one, and output read list; 
# Step 3 : Retrieve read alignments from transRdBam; 
use strict; 
use warnings; 
use LogInforSunhh;
use SeqAlnSunhh;
use fileSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts, 
	"help!", 
	"pref:s",    # Project prefix; 
	"wrk_dir:s", # './', Assign a directory to store all resulting bam files, instead of in current folder.
	"transRdBam:s@", # Input pref.src2tgt_rd.bam ; 
	"bgRdBam:s@",    # Input pref_comb.nonTransRmdup.bam ; 
	"tgt_fa:s@",     # Input target_genome.ref.fasta ; 
	"outBam:s",      # Output filtered bam; pref.src2tgt_cleanRd.bam
	"transRdMaxLen:i", # cut transmitted read to shorter length to fit -bgRdBam 

	"exe_samtools:s", 
	"pl_getAln:s",   # $gg{'dir_abs'}/get_alnBam_by_src2tgt_rdList.pl 

	"startStep:i", 
); 

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) };

my %gg; 
$gg{'windSize'} = 1e6; 
$gg{'isSrcRdID'}   = {}; # Read here has been checked. 1-src, 0-tgt; 
$gg{'isSrcSeq'}   = {}; # Read sequence here has been checked. 1-src, 0-tgt; 

&setGlob(); 
&applyOpt(); 
&step3_cleanByTgtRd(); 


sub step3_cleanByTgtRd {
	# Step 3 : Clean possbile transmitted reads by reads from target-only samples. 
	# Produce : pref.src2tgt_cleanRd.bam : $gg{'wrk_dir'}/$gg{'pref'}.src2tgt_cleanRd.bam
	
	$gg{'startStep'} <= 0 and do { scalar(@{$gg{'transRdBam'}}) > 0 or &stopErr("[Err] At least one -transRdBam is required.\n") }; 
	for my $fn (@{$gg{'transRdBam'}}) {
		$gg{'startStep'} <= 0 and do { -e $fn or &stopErr("[Err] Failed to find file -transRdBam [$fn]\n"); }; 
	}
	my @toRM; 
	
	# Format bam files. 
	my $bamlist = "$gg{'wrk_dir'}/$gg{'pref'}.step3.iBamList"; 
	my $bamlistSrc = "$gg{'wrk_dir'}/$gg{'pref'}.step3.iBamListSrc"; 
	my $cmd = ""; 
	$gg{'startStep'} <= 0 and &fileSunhh::write2file($bamlist,'','>'); 
	$gg{'startStep'} <= 0 and &fileSunhh::write2file($bamlistSrc, '', '>'); 
	for (my $i=0; $i<@{$gg{'transRdBam'}}; $i++) {
		my $newFn = "$gg{'wrk_dir'}/$gg{'pref'}.src_$i.bam"; 
		$gg{'startStep'} <= 0 and &fmtBamRG( $gg{'transRdBam'}[$i], $newFn, "src" ); 
		$gg{'startStep'} <= 0 and &fileSunhh::write2file($bamlist, "$newFn\n",'>>'); 
		$gg{'startStep'} <= 0 and &fileSunhh::write2file($bamlistSrc, "$gg{'transRdBam'}[$i]", '>>'); 
		$gg{'startStep'} <= 0 and push(@toRM, $newFn); 
	}
	for (my $i=0; $i<@{$gg{'bgRdBam'}}; $i++) {
		my $newFn = "$gg{'wrk_dir'}/$gg{'pref'}.tgt_$i.bam"; 
		$gg{'startStep'} <= 0 and &fmtBamRG( $gg{'bgRdBam'}[$i], $newFn, "tgt" ); 
		$gg{'startStep'} <= 0 and &fileSunhh::write2file($bamlist, "$newFn\n",'>>'); 
		$gg{'startStep'} <= 0 and push(@toRM, $newFn); 
	}
	$gg{'startStep'} <= 0 and push(@toRM, $bamlist, $bamlistSrc); 
	# Combine bam files. 
	$gg{'startStep'} <= 0 and $cmd .= "$gg{'exe_samtools'} merge -c -f "; 
	$gg{'startStep'} <= 0 and $cmd .= " -\@ 10 "; 
	$gg{'startStep'} <= 0 and $cmd .= " -b $bamlist "; 
	$gg{'startStep'} <= 0 and $cmd .= " $gg{'wrk_dir'}/$gg{'pref'}.step3.jn.bam "; 
	$gg{'startStep'} <= 0 and &exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$gg{'startStep'} <= 0 and $cmd = ""; 
	$gg{'startStep'} <= 0 and push(@toRM, "$gg{'wrk_dir'}/$gg{'pref'}.step3.jn.bam"); 

	# Sort joined bam file. 
	$gg{'startStep'} <= 1 and $cmd .= "$gg{'exe_samtools'} sort "; 
	$gg{'startStep'} <= 1 and $cmd .= " -\@ 10 -m 5G "; 
	$gg{'startStep'} <= 1 and $cmd .= " -o $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam "; 
	$gg{'startStep'} <= 1 and $cmd .= " $gg{'wrk_dir'}/$gg{'pref'}.step3.jn.bam "; 
	$gg{'startStep'} <= 1 and &exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$gg{'startStep'} <= 1 and $cmd = ""; 
	$gg{'startStep'} <= 1 and push(@toRM, "$gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam"); 

	# Check reads by sliding window; 
	my %curH; 
	my @notTgtRd; 
if ( $gg{'startStep'} <= 2 ) {
	open F,'-|', "$gg{'exe_samtools'} view $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam" or &stopErr("[Err] Failed to read bam file [$gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam]\n"); 
	my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>1e5); 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		defined $flag_UN{$ta[1]} and next; 
		$ta[5] eq '*' and next; 
		&fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Processing $. line for $gg{'pref'}.\n"); 
		if (defined $curH{'refID'}) {
			if ($curH{'refID'} ne $ta[2]) {
				my ($notTgtRd_ar) = &process_curH(\%curH, 'end'); 
				push(@notTgtRd, @$notTgtRd_ar); 
				&add2curH(\%curH, \@ta, 'y'); 
			} elsif ( $ta[3] > $curH{'startP'} + $gg{'windSize'} - 1 ) {
				my ($notTgtRd_ar) = &process_curH(\%curH); 
				push(@notTgtRd, @$notTgtRd_ar); 
				&add2curH(\%curH, \@ta, 'n'); 
			} else {
				&add2curH(\%curH, \@ta, 'n'); 
			}
		} else {
			&add2curH(\%curH, \@ta, 'y'); 
		}
	}
	close F; 
	my ($notTgtRd_ar) = &process_curH(\%curH, 'end'); 
	push(@notTgtRd, @$notTgtRd_ar); 
	my $ofhRd = &openFH("$gg{'wrk_dir'}/$gg{'pref'}.step3.notTgtRd.list", '>'); 
	my %hasOutRdID; 
	for (@notTgtRd) {
		$gg{'isSrcRdID'}{$_} == 1 or next; 
		defined $hasOutRdID{$_} and next; 
		$hasOutRdID{$_} = 1; 
		print {$ofhRd} "$_\n"; 
	}
	close ($ofhRd); 
}# End if ($gg{'startStep'} <= 2)
	$gg{'startStep'} <= 3 and $cmd = ""; 
	$gg{'startStep'} <= 3 and $cmd .= "$gg{'exe_samtools'} merge -c -f "; 
	$gg{'startStep'} <= 3 and $cmd .= " -\@ 10 "; 
	$gg{'startStep'} <= 3 and $cmd .= " -b $bamlistSrc "; 
	$gg{'startStep'} <= 3 and $cmd .= " $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrc.bam "; 
	$gg{'startStep'} <= 3 and &exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$gg{'startStep'} <= 3 and $cmd = ""; 
	$gg{'startStep'} <= 3 and push(@toRM, "$gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrc.bam"); # Commended temp; 
	$gg{'startStep'} <= 4 and $cmd .= "perl $gg{'pl_getAln'} "; 
	$gg{'startStep'} <= 4 and $cmd .= " -wrk_dir  $gg{'wrk_dir'} -pref $gg{'pref'}.step3 "; 
	$gg{'startStep'} <= 4 and $cmd .= " -inBam    $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrc.bam "; 
	$gg{'startStep'} <= 4 and $cmd .= " -inRdList $gg{'wrk_dir'}/$gg{'pref'}.step3.notTgtRd.list "; 
	$gg{'startStep'} <= 4 and $cmd .= " -outBam   $gg{'outBam'}  "; 
	$gg{'startStep'} <= 4 and $cmd .= " -exe_samtools $gg{'exe_samtools'} "; 
	$gg{'startStep'} <= 4 and &exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$gg{'startStep'} <= 4 and $cmd = ""; 

	for (@toRM) {
		$gg{'startStep'} <= 4 and &fileSunhh::_rmtree($_); 
	}

	return; 
}# step3_cleanByTgtRd() 
sub process_curH {
	my ($cH, $all) = @_; 
	$all //= 'n'; # If $all eq 'n', this is not the end of the chromosome, or else this %curH has all the alignments at the end of this CHR. 

	my $max_span = 200e3; 
	my (@notTgtRd); 
	my @allAln = @{$cH->{'aln'}}; 
	my $cntAln = scalar(@allAln); 
	my $stopChkS = $cH->{'startP'} + $gg{'windSize'} - $max_span; 
	$stopChkS < $cH->{'startP'} and &stopErr("[Err] StopChkS [$stopChkS] smaller than start [$cH->{'startP'}]\n"); 
	$all eq 'n' or $stopChkS = $allAln[-1][2]+$max_span; 
	&tsmsg("[Msg] Processing $cH->{'refID'} $cH->{'startP'} to $stopChkS with $cntAln alignments for $gg{'pref'}\n"); 
	$cntAln > 0 or return(\@notTgtRd); 
	my @tgtAlnIdx; 
	my %srcAlnIdx2tgtIdx; 
	for (my $i=0; $i<@allAln; $i++) {
		if ( $allAln[$i][3] == 1 ) {
			# This read comes from src; 
			if ( $#tgtAlnIdx > -1 ) {
				$srcAlnIdx2tgtIdx{$i} = [$#tgtAlnIdx, $#tgtAlnIdx+1]; 
			} else {
				$srcAlnIdx2tgtIdx{$i} = [-1, 0]; 
			}
		} else {
			push(@tgtAlnIdx, $i); 
		}
	}
	&tsmsg("[Msg]   There are " . scalar(@tgtAlnIdx) . " TGT-alignments for $gg{'pref'}\n"); 
	&tsmsg("[Msg]   There are " . scalar(keys %srcAlnIdx2tgtIdx) . " SRC-alignments for $gg{'pref'}\n"); 

	my (@newAln); 
	my $i=-1; 
	for my $cH_curAln (@allAln) {
		$i++; 
		$i % 1000 == 0 and &tsmsg("[Msg]   Processing $i -th alignment for $gg{'pref'} within [$cH->{'refID'} $cH->{'startP'} to $stopChkS]\n"); 
		if ( $cH_curAln->[1] > $stopChkS ) {
			push(@newAln, $cH_curAln); 
			next; 
		}
		# $cH_curAln->[1] >= $beginNew and push(@newAln, $cH_curAln); 
		# $cH_curAln->[1] > $stopChkS or next; 
		$cH_curAln->[3] == 1 or next; # Only SRC read is considered. 
		my $tk_h_seq = "$cH_curAln->[1]\t$cH_curAln->[4]"; # "startPos\trdSeq\n"; 
		my $tk_h_subseq = $tk_h_seq; 
		if ( defined $gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq} ) {
			if ( $gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq} == 1 ) {
				# Value '1' means this sequence mapping here is a real src2tgt read. 
				if (defined $gg{'isSrcRdID'}{ $cH_curAln->[0] }) {
					if ( $gg{'isSrcRdID'}{ $cH_curAln->[0] } == 1 ) {
						; 
					} else {
						# $gg{'isSrcRdID'}{ $cH_curAln->[0] } == 0 
						# I should correct the bad classification of $gg{'hasChkSeq'}{$cH->{'refID'}}{$tk_h_seq}; 
						$gg{'hasChkSeq'}{$cH->{'refID'}}{$tk_h_seq} = 0; 
					}
				} else {
					# This read has not been checked. 
					push(@notTgtRd, $cH_curAln->[0]); 
					$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 1; 
				}
			} else {
				# $gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq} == 0, should be included by TGT read. 
				$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 0; 
			}
			next; 
		}
		if ( defined $gg{'isSrcRdID'}{ $cH_curAln->[0] } and $gg{'isSrcRdID'}{ $cH_curAln->[0] } == 0 ) {
			# Should be included by TGT read. 
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq} = 0; 
			next; 
		}
		# This read might be 'isSrcRdID'==1, or not checked yet. 
		# I am thinking if a read is 'isSrcRdID'==1, is it possible it could be changed to 'isSrcRdID'=0 by another alignment. 
		# defined $gg{'isSrcRdID'}{ $cH_curAln->[0] } and next; 
		my $is_tgt = 0; 
		my $subseq = $cH_curAln->[4]; 
		my $subseqL = $cH_curAln->[5]; 
		if ( $gg{'transRdMaxLen'} > 1 and $gg{'transRdMaxLen'} < $subseqL ) {
			$subseq  = substr( $subseq, 0, $gg{'transRdMaxLen'} ); 
			$subseqL = length( $subseq ); 
			$tk_h_subseq = "$cH_curAln->[1]\t$subseq"; # "startPos\trdSeq\n"; 
			if ( defined $gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} ) {
				if ( $gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} == 1 ) {
					# Value '1' means this sequence mapping here is a real src2tgt read. 
					if (defined $gg{'isSrcRdID'}{ $cH_curAln->[0] }) {
						if ( $gg{'isSrcRdID'}{ $cH_curAln->[0] } == 1 ) {
							; 
						} else {
							# $gg{'isSrcRdID'}{ $cH_curAln->[0] } == 0 
							# I should correct the bad classification of $gg{'hasChkSeq'}{$cH->{'refID'}}{$tk_h_subseq}; 
							$gg{'hasChkSeq'}{$cH->{'refID'}}{$tk_h_subseq} = 0; 
						}
					} else {
						# This read has not been checked. 
						push(@notTgtRd, $cH_curAln->[0]); 
						$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 1; 
					}
				} else {
					# $gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} == 0, should be included by TGT read. 
					$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 0; 
				}
				next; 
			}
			if ( defined $gg{'isSrcRdID'}{ $cH_curAln->[0] } and $gg{'isSrcRdID'}{ $cH_curAln->[0] } == 0 ) {
				# Should be included by TGT read. 
				$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} = 0; 
				next; 
			}
		}

		for (my $j0=$srcAlnIdx2tgtIdx{$i}[0]; $j0>=0; $j0--) {
			my $t = $allAln[ $tgtAlnIdx[$j0] ]; 
			# $t->[3] == 0 or next; # Only 0-tgt is used to check TGT; 
			$t->[1] + $max_span < $cH_curAln->[1] and last; 
			$t->[2] < $cH_curAln->[1] and next; 
			$t->[5] < $subseqL and next; 
			index($t->[4], $subseq) != -1 and do { $is_tgt = 1; last; }; 
		}
		if ($is_tgt == 1) {
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq}    = 0; # Value '0' means this sequence mapping here is included by a TGT read. 
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} = 0; 
			$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 0; 
			next; 
		}
		for (my $j0=$srcAlnIdx2tgtIdx{$i}[1]; $j0<@tgtAlnIdx; $j0++) {
			my $t = $allAln[ $tgtAlnIdx[$j0] ]; 
			$t->[3] == 0 or next; # Only 0-tgt is used to check TGT; 
			$t->[1] > $cH_curAln->[1] and last; 
			$t->[5] < $subseqL and next; 
			index($t->[4], $subseq) != -1 and do { $is_tgt = 1; last; }; 
		}
		if ($is_tgt == 1) {
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq}    = 0; # Value '0' means this sequence mapping here is included by a TGT read. 
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} = 0; 
			$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 0; 
			next; 
		}
		if ( defined $gg{'isSrcRdID'}{ $cH_curAln->[0] } ) {
			# This read has been checked, and classified as 'isSrcRdID'==1; 
			next; 
		} else {
			push(@notTgtRd, $cH_curAln->[0]); 
			$gg{'isSrcRdID'}{ $cH_curAln->[0] } = 1; 
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_seq}    = 1; # Value '1' means this sequence mapping here is not included by a TGT read. 
			$gg{'isSrcSeq'}{$cH->{'refID'}}{$tk_h_subseq} = 1; 
		}
	}

	@{$cH->{'aln'}} = @newAln; 
	$cH->{'startP'} = $stopChkS+1; 
	my @back_notTgtRd; 
	for my $t (@notTgtRd) {
		$gg{'isSrcRdID'}{$t} == 1 and push(@back_notTgtRd, $t); 
	}

	return(\@back_notTgtRd); 
}# process_curH() 

sub add2curH {
	my ($cH, $taR, $first) = @_; 
	$first //= 'n'; 
	if ($first ne 'n') {
		$cH->{'refID'}  = $taR->[2]; 
		$cH->{'startP'} = $taR->[3]; 
	}
	# $taR->[11] eq 'RG:Z:src' and $cH->{'srcN'} ++; 
	my ($rdLen, $ref_span) = &SeqAlnSunhh::cigar_array2len( &SeqAlnSunhh::cigar_str2array( $taR->[5] ) ); 
	my @o = ('', '', '', '', '', ''); 
	# 0          1          2                      3            4          5
	# rdID     , startPos , endPos               , Src_1/Tgt_0, rdSeq    , rdLength
	if      ( $taR->[11] eq 'RG:Z:tgt' ) {
		@o = ('',        $taR->[3], $taR->[3]+$ref_span-1, 0,           $taR->[9], length($taR->[9])); 
	} elsif ( $taR->[11] eq 'RG:Z:src' ) {
		@o = ($taR->[0], $taR->[3], $taR->[3]+$ref_span-1, 1,           $taR->[9], length($taR->[9])); 
		#     0          1          2                      3            4          5
		#     rdID     , startPos , endPos               , Src_1/Tgt_0, rdSeq    , rdLength
	} else {
		&stopErr("[Err] Bad input aln line with [$taR->[11]]: @$taR\n"); 
	}
	push(@{$cH->{'aln'}}, \@o); 
	return; 
}# add2curH()

sub add2curH_0 {
	# Deprecated. 
	my ($cH, $taR, $first) = @_; 
	$first //= 'n'; 
	if ($first ne 'n') {
		$cH->{'refID'}  = $taR->[2]; 
		$cH->{'startP'} = $taR->[3]; 
	}
	# $taR->[11] eq 'RG:Z:src' and $cH->{'srcN'} ++; 
	my ($rdLen, $ref_span) = &SeqAlnSunhh::cigar_array2len( &SeqAlnSunhh::cigar_str2array( $taR->[5] ) ); 
	push(@{$cH->{'aln'}}, [@$taR, $taR->[3]+$ref_span-1]); 
}# add2curH_0()

sub fmtBamRG {
	my ($iBam, $oBam, $rgID) = @_; 
	&tsmsg("[Msg] Reformatting bam file [$iBam]\n"); 
	open F,'-|', "$gg{'exe_samtools'} view -h $iBam " or &stopErr("[Err] Failed to read bam file [$iBam]\n");
	open O,'|-', "$gg{'exe_samtools'} view -o $oBam - " or &stopErr("[Err] Failed to write bam file [$oBam]\n"); 
	my $has_RG = 0; 
	while (<F>) {
		if (m!^\@!) {
			m!^\@(RG|PG)\t! and next; 
			print O $_; 
			next; 
		}
		$has_RG == 0 and do { $has_RG = 1; print O join("\t", '@RG', "ID:$rgID", "SM:$rgID", "LB:$rgID", "PL:illumina")."\n"; }; 
		m!^([^\t]+(?:\t[^\t]+){10})(\t|$)! or &stopErr("[Err] Bad format line: $_\n"); 
		print O "$1\tRG:Z:$rgID\n"; 
	}
	close O; 
	close F; 

	return; 
}#fmtBamRG() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 

	# For outer tools
	for my $k (qw/exe_samtools pl_getAln/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	for my $k (qw/pref wrk_dir transRdBam bgRdBam tgt_fa outBam transRdMaxLen startStep/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	defined $opts{'transRdBam'} or $gg{'transRdBam'} = [ "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rd.bam" ]; 
	defined $opts{'bgRdBam'} or $gg{'bgRdBam'} = [ "$gg{'wrk_dir'}/$gg{'pref'}.nonTransRmdup.bam" ]; 
	defined $opts{'tgt_fa'} or $gg{'tgt_fa'} = [ ]; 
	defined $opts{'outBam'} or $gg{'outBam'} = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_cleanRd.bam"; 

	return; 
}# applyOpt() 

sub setGlob {
	$gg{'pref'}    = 'pref'; 
	$gg{'wrk_dir'} = "./"; 

	$gg{'transRdBam'} = [ "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_rd.bam" ]; 
	$gg{'bgRdBam'}    = [ "$gg{'wrk_dir'}/$gg{'pref'}.nonTransRmdup.bam" ]; 
	$gg{'tgt_fa'}     = [ '' ]; 
	$gg{'outBam'}  = "$gg{'wrk_dir'}/$gg{'pref'}.src2tgt_cleanRd.bam"; 

	$gg{'exe_samtools'} = 'samtools'; 

	$gg{'abs_cur_dir'} = &fileSunhh::_abs_path("./"); 
	$gg{'abs_wrk_dir'} = &fileSunhh::_abs_path($gg{'wrk_dir'}); 
	$gg{'abs_dir_basePL'}  = &fileSunhh::_dirname(&fileSunhh::_abs_path($0)); 
	$gg{'link_dir_basePL'} = &fileSunhh::_dirname(&fileSunhh::_abs_path_4link($0)); 
	$gg{'transRdMaxLen'}   = -1; 

	$gg{'pl_getAln'}  = ( -e "$gg{'abs_dir_basePL'}/get_alnBam_by_src2tgt_rdList.pl" ) ? "$gg{'abs_dir_basePL'}/get_alnBam_by_src2tgt_rdList.pl" : "$gg{'link_dir_basePL'}/get_alnBam_by_src2tgt_rdList.pl" ; 

	$gg{'startStep'}  = 0; 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# step 3 . Check if possible transmitted reads are contained by non/self-grafted samples.
#
# perl $0   -pref $gg{'pref'} -wrk_dir $gg{'wrk_dir'}
#
#   -pref         [$gg{'pref'}] Output prefix ; 
#     -outBam     [$gg{'outBam'}] Output filtered bam; pref.src2tgt_cleanRd.bam
#     -transRdBam [$gg{'transRdBam'}[0]] \@ Input pref.src2tgt_rd.bam ;
#   -transRdMaxLen [$gg{'transRdMaxLen'}] cut transmitted read to shorter length to fit -bgRdBam
#   -bgRdBam      [$gg{'bgRdBam'}[0]] \@ Input pref_comb.nonTransRmdup.bam ;
#     -tgt_fa     [$gg{'tgt_fa'}[0]] \@ Not necessary. Input target_genome.ref.fasta ; 
#   -wrk_dir      [$gg{'wrk_dir'}] Directory storing result files. 
#
#   -exe_samtools      [$gg{'exe_samtools'}]
#   -pl_getAln         [$gg{'pl_getAln'}]
################################################################################
HH
	return; 
}# setGlob() 


