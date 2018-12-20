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
	
	"testChk!", 
); 

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) };

my %gg; 
$gg{'windSize'} = 1e6; 
&setGlob(); 
&applyOpt(); 
&step3_cleanByTgtRd(); 


sub step3_cleanByTgtRd {
	# Step 3 : Clean possbile transmitted reads by reads from target-only samples. 
	# Produce : pref.src2tgt_cleanRd.bam : $gg{'wrk_dir'}/$gg{'pref'}.src2tgt_cleanRd.bam
	
	scalar(@{$gg{'transRdBam'}}) > 0 or &stopErr("[Err] At least one -transRdBam is required.\n"); 
	for my $fn (@{$gg{'transRdBam'}}) {
		-e $fn or &stopErr("[Err] Failed to find file -transRdBam [$fn]\n"); 
	}
	my @toRM; 
	
	# Format bam files. 
	my $bamlist = "$gg{'wrk_dir'}/$gg{'pref'}.step3.iBamList"; 
	my $bamlistSrc = "$gg{'wrk_dir'}/$gg{'pref'}.step3.iBamListSrc"; 
	my $cmd = ""; 
if (!$opts{'testChk'}) {
	&fileSunhh::write2file($bamlist,'','>'); 
	&fileSunhh::write2file($bamlistSrc, '', '>'); 
	for (my $i=0; $i<@{$gg{'transRdBam'}}; $i++) {
		my $newFn = "$gg{'wrk_dir'}/$gg{'pref'}.src_$i.bam"; 
		&fmtBamRG( $gg{'transRdBam'}[$i], $newFn, "src" ); 
		&fileSunhh::write2file($bamlist, "$newFn\n",'>>'); 
		&fileSunhh::write2file($bamlistSrc, "$gg{'transRdBam'}[$i]", '>>'); 
		push(@toRM, $newFn); 
	}
	for (my $i=0; $i<@{$gg{'bgRdBam'}}; $i++) {
		my $newFn = "$gg{'wrk_dir'}/$gg{'pref'}.tgt_$i.bam"; 
		&fmtBamRG( $gg{'bgRdBam'}[$i], $newFn, "tgt" ); 
		&fileSunhh::write2file($bamlist, "$newFn\n",'>>'); 
		push(@toRM, $newFn); 
	}
	push(@toRM, $bamlist, $bamlistSrc); 
	# Combine bam files. 
	$cmd .= "$gg{'exe_samtools'} merge -c -f "; 
	$cmd .= " -\@ 10 "; 
	$cmd .= " -b $bamlist "; 
	$cmd .= " $gg{'wrk_dir'}/$gg{'pref'}.step3.jn.bam "; 
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$cmd = ""; 
	push(@toRM, "$gg{'wrk_dir'}/$gg{'pref'}.step3.jn.bam"); 

	# Sort joined bam file. 
	$cmd .= "$gg{'exe_samtools'} sort "; 
	$cmd .= " -\@ 10 -m 5G "; 
	$cmd .= " -o $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam "; 
	$cmd .= " $gg{'wrk_dir'}/$gg{'pref'}.step3.jn.bam "; 
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$cmd = ""; 
	push(@toRM, "$gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam"); 

}
	# Check reads by sliding window; 
	my %curH; 
	my @notTgtRd; 
	open F,'-|', "$gg{'exe_samtools'} view $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam" or &stopErr("[Err] Failed to read bam file [$gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrt.bam]\n"); 
	my %tmp_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>1e5); 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		defined $flag_UN{$ta[1]} and next; 
		$ta[5] eq '*' and next; 
		&fileSunhh::log_section( $. , \%tmp_cnt ) and &tsmsg("[Msg] Processing $. line.\n"); 
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
	for (@notTgtRd) {
		print {$ofhRd} "$_\n"; 
	}
	close ($ofhRd); 
	$cmd = ""; 
	$cmd .= "$gg{'exe_samtools'} merge -c -f "; 
	$cmd .= " -\@ 10 "; 
	$cmd .= " -b $bamlistSrc "; 
	$cmd .= " $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrc.bam "; 
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$cmd = ""; 
	push(@toRM, "$gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrc.bam"); 
	$cmd .= "perl $gg{'pl_getAln'} "; 
	$cmd .= " -wrk_dir  $gg{'wrk_dir'} -pref $gg{'pref'}.step3 "; 
	$cmd .= " -inBam    $gg{'wrk_dir'}/$gg{'pref'}.step3.jnSrc.bam "; 
	$cmd .= " -inRdList $gg{'wrk_dir'}/$gg{'pref'}.step3.notTgtRd.list "; 
	# $cmd .= " -outBam   $gg{'wrk_dir'}/$gg{'pref'}.step3.notTgtRd.bam "; 
	$cmd .= " -outBam   $gg{'outBam'}  "; 
	$cmd .= " -exe_samtools $gg{'exe_samtools'} "; 
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	$cmd = ""; 

	for (@toRM) {
		&fileSunhh::_rmtree($_); 
	}

	return; 
}# step3_cleanByTgtRd() 
sub process_curH {
	my ($cH, $all) = @_; 
	$all //= 'n'; # If $all eq 'n', this is not the end of the chromosome, or else this %curH has all the alignments at the end of this CHR. 

	my (@notTgtRd, %h); 
	scalar($cH->{'aln'}) > 0 or return(\@notTgtRd); 
	my $max_span = 200e3; 

	my $stopChkS = $cH->{'startP'} + $gg{'windSize'} - $max_span; 
	$stopChkS < $cH->{'startP'} and &stopErr("[Err] StopChkS [$stopChkS] smaller than start [$cH->{'startP'}]\n"); 
	$all eq 'n' or $stopChkS = $cH->{'aln'}[-1][3]+$max_span; 
	# $opts{'testChk'} and &tsmsg("[Msg] Processing $cH->{'refID'} $cH->{'startP'} to $stopChkS\n"); 
	&tsmsg("[Msg] Processing $cH->{'refID'} $cH->{'startP'} to $stopChkS\n"); 
	my (@newAln); 
	for (my $i=0; $i<@{$cH->{'aln'}}; $i++) {
		if ( $cH->{'aln'}[$i][3] > $stopChkS ) {
			push(@newAln, [ @{$cH->{'aln'}[$i]} ]); 
			next; 
		}
		# $cH->{'aln'}[$i][3] >= $beginNew and push(@newAln, [ @{$cH->{'aln'}[$i]} ]); 
		# $cH->{'aln'}[$i][3] > $stopChkS or next; 
		$cH->{'aln'}[$i][11] eq 'RG:Z:src' or next; 
		defined $cH->{'chked'}{ $cH->{'aln'}[$i][0] } and next; 
		$cH->{'chked'}{ $cH->{'aln'}[$i][0] } = 1; 
		my $is_tgt = 0; 
		my $subseq = $cH->{'aln'}[$i][9]; 
		$gg{'transRdMaxLen'} > 1 and length($subseq) > $gg{'transRdMaxLen'} and $subseq = substr( $subseq, 0, $gg{'transRdMaxLen'} ); 
		for (my $j=$i-1; $j>=0; $j--) {
			$cH->{'aln'}[$j][11] eq 'RG:Z:src' and next; 
			$cH->{'aln'}[$j][3]+$max_span < $cH->{'aln'}[$i][3] and last; 
			$cH->{'aln'}[$j][12] < $cH->{'aln'}[$i][3] and next; 
			index($cH->{'aln'}[$j][9], $subseq) != -1 and do { $is_tgt = 1; last; }; 
		}
		$is_tgt == 1 and next; 
		for (my $j=$i+1; $j<@{$cH->{'aln'}}; $j++) {
			$cH->{'aln'}[$j][11] eq 'RG:Z:src' and next; 
			$cH->{'aln'}[$j][3] > $cH->{'aln'}[$i][3] and last; 
			index($cH->{'aln'}[$j][9], $subseq) != -1 and do { $is_tgt = 1; last; }; 
		}
		$is_tgt == 1 and next; 
		defined $h{ $cH->{'aln'}[$i][0] } or push(@notTgtRd, $cH->{'aln'}[$i][0]); 
		$cH->{'cntNotTgtRd'} ++; 
		$h{ $cH->{'aln'}[$i][0] } = 1; 
	}

	@{$cH->{'aln'}} = @newAln; 
	$cH->{'startP'} = $stopChkS+1; 

	return(\@notTgtRd); 
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
	push(@{$cH->{'aln'}}, [@$taR, $taR->[3]+$ref_span-1]); 
}# add2curH()

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

	for my $k (qw/pref wrk_dir transRdBam bgRdBam tgt_fa outBam transRdMaxLen/) {
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


