#!/usr/bin/perl
# 0 : unmapped ; 1 : mapped with mismatch ; 2 : mapped with 0-mismatch 
# Separate input reads into XX groups : 
#   Group 1 : map with 0-mismatch in both SAMs; 
#   Group 2 : map with 0-mismatch in SAM1 and with mismatch in SAM2 ; 
#   Group 3 : map with 0-mismatch in SAM2 and with mismatch in SAM1 ; 
#   Group 4 : map with 0-mismatch in SAM1 and not mapped in SAM2 ; 
#   Group 5 : map with 0-mismatch in SAM2 and not mapped in SAM1 ; 
#   Group 6 : map with mismatch in both SAMs ; 
#   Group 7 : map with mismatch in SAM1 and not mapped in SAM2 ; 
#   Group 8 : map with mismatch in SAM2 and not mapped in SAM1 ; 
#   Group 9 : Others including - 'not mapped in any SAMs' ; 
# If input reads are paired, two ends are considered together. 
#   Mapping : both ends can be mapped ; 
#   0-mismatch : both ends are mapped with 0-mismatch ; 
# For mismatch calculation : 
#   Cigars [HSIDX] are all considered as mismatches. 
#   Because I use hisat2, so 'XM:i:n' is also considered. 
#
use strict; 
use warnings; 
use fileSunhh; 
use mathSunhh; 
use SeqAlnSunhh; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"tell_grp!", 
	"verbose!", 
	"inSam1:s", "inSam2:s", 
	"inFq1:s",  "inFq2:s", 
	"fmtSam1:s", # sam
	"fmtSam2:s", # sam
	"outPref:s", # out

	"exe_samtools:s", # samtools 
	"exe_perl:s", # perl 
	"pl_extractFq:s", # /home/Sunhh/tools/github/NGS_data_processing/extract_fq_by_list.pl

); 

my %grpSE; 
%grpSE = qw(
2_2   1
2_1   2
1_2   3
2_0   4
0_2   5
1_1   6
1_0   7
0_1   8
0_0   9
); 

my %grpPE; 
%grpPE = qw(
2_2   1
2_1   2
1_2   3
2_0   4
0_2   5
1_1   6
1_0   7
0_1   8
0_0   9
); 



my $help_txt = <<HH; 

perl $0 -inSam1 in_toRef1.sam   -inSam2 in_toRef2.sam   -inFq1 in_src.fq [ -inFq2 in_src_R2.fq ]

-outPref      ['out']

-fmtSam1      ['sam' or detected]
-fmtSam2      ['sam' or detected]

-exe_samtools ['samtools']
-exe_perl     ['perl']
-pl_extractFq [/home/Sunhh/tools/github/NGS_data_processing/extract_fq_by_list.pl]

-help         [Boolean]
-tell_grp     [Boolean]

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
$opts{'tell_grp'} and &tell_info(); 


my %cnt; 

&prepare_input(); 

my %rd_info = %{ &load_fqID( $opts{'inFq1'} ) }; 
#   return({ 'ID2idx'=>\%back_hash, 'map_idx'=>[@back_id] });
$cnt{'inRdNum'} = scalar(@{$rd_info{'map_idx'}}); 
$opts{'verbose'} and &tsmsg("[Msg] inFq1 reads number : $cnt{'inRdNum'}\n"); 
&set_info_bySam( $opts{'inSam1'}, $opts{'fmtSam1'}, 1, \%rd_info ); 
&set_info_bySam( $opts{'inSam2'}, $opts{'fmtSam2'}, 2, \%rd_info ); 

my $wrk_dir = &fileSunhh::new_tmp_dir(); 
mkdir($wrk_dir) or &stopErr("[Err] Failed to create working directory [$wrk_dir]\n"); 
if ( defined $opts{'inFq2'} ) {
	$opts{'verbose'} and &tsmsg("[Msg] Grouping reads.\n"); 
	for my $ar1 (@{ $rd_info{'map_idx'} }) {
		my $min1 = &mathSunhh::min( $ar1->[0], $ar1->[1] ); 
		my $min2 = &mathSunhh::min( $ar1->[2], $ar1->[3] ); 
		my $tag = "${min1}_${min2}"; 
		my $grpID = ( defined $grpPE{$tag} ) ? $grpPE{$tag} : 9 ; 
		&fileSunhh::write2file( "$wrk_dir/grp$grpID.rdID", "$ar1->[4]\n", '>>' ); 
		$cnt{'RdNInGrp'}{$grpID} ++; 
	}
	for my $grpID (sort {$a<=>$b} keys %{ { map { $_ => 1 } values %grpSE } }) {
		$opts{'verbose'} and &tsmsg("[Msg] Extracting fastq for grp[$grpID]\n"); 
		$cnt{'RdNInGrp'}{$grpID} //= 0; 
		&tsmsg( join('', "[Msg] Grp${grpID} pairs number : $cnt{'RdNInGrp'}{$grpID} (", sprintf("%0.4f", $cnt{'RdNInGrp'}{$grpID}/$cnt{'inRdNum'}*100),"%)\n") ); 
		if ( $cnt{'RdNInGrp'}{$grpID} == 0 ) {
			&tsmsg("[Msg] No reads in grp$grpID\n"); 
			&exeCmd_1cmd("rm -f $opts{'outPref'}R1.grp$grpID.fq ; touch $opts{'outPref'}R1.grp$grpID.fq"); 
			&exeCmd_1cmd("rm -f $opts{'outPref'}R2.grp$grpID.fq ; touch $opts{'outPref'}R2.grp$grpID.fq"); 
		} else {
			&exeCmd_1cmd( "$opts{'exe_perl'} $opts{'pl_extractFq'} -mode keep -rdKey -trim12   -refLis $wrk_dir/grp$grpID.rdID   -srcFq $opts{'inFq1'}   -outFq $opts{'outPref'}R1.grp$grpID.fq" ); 
			&exeCmd_1cmd( "$opts{'exe_perl'} $opts{'pl_extractFq'} -mode keep -rdKey -trim12   -refLis $wrk_dir/grp$grpID.rdID   -srcFq $opts{'inFq2'}   -outFq $opts{'outPref'}R2.grp$grpID.fq" ); 
		}
	}
} else {
	$opts{'verbose'} and &tsmsg("[Msg] Grouping reads.\n"); 
	for my $ar1 (@{ $rd_info{'map_idx'} }) {
		my $tag = "$ar1->[0]_$ar1->[2]"; 
		# my $grpID = ( defined $grpSE{$tag} ) ? $grpSE{$tag} : 9 ; 
		my $grpID = $grpSE{$tag}; 
		&fileSunhh::write2file( "$wrk_dir/grp$grpID.rdID", "$ar1->[4]\n", '>>' ); 
		$cnt{'RdNInGrp'}{$grpID} ++; 
	}
	for my $grpID (sort {$a<=>$b} keys %{ { map { $_ => 1 } values %grpSE } }) {
		$cnt{'RdNInGrp'}{$grpID} //= 0; 
		&tsmsg( join('', "[Msg] Grp${grpID} pairs number : $cnt{'RdNInGrp'}{$grpID} (", sprintf("%0.4f", $cnt{'RdNInGrp'}{$grpID}/$cnt{'inRdNum'}*100),"%)\n") ); 
		$opts{'verbose'} and &tsmsg("[Msg] Extracting fastq for grp[$grpID]\n"); 
		if ( $cnt{'RdNInGrp'}{$grpID} == 0 ) {
			&tsmsg("[Msg] No reads in grp$grpID\n"); 
			&exeCmd_1cmd("rm -f $opts{'outPref'}grp$grpID.fq ; touch $opts{'outPref'}grp$grpID.fq"); 
		} else {
			&exeCmd_1cmd( "$opts{'exe_perl'} $opts{'pl_extractFq'} -mode keep -rdKey -trim12   -refLis $wrk_dir/grp$grpID.rdID   -srcFq $opts{'inFq1'}   -outFq $opts{'outPref'}grp$grpID.fq" ); 
		}
	}
}
&fileSunhh::_rmtree( $wrk_dir ); 

&tsmsg("[Rec] script [$0] done.\n"); 

sub set_info_bySam {
	my ( $fnSam, $fmtSam, $samOrder, $rd_info_h ) = @_; 
	$opts{'verbose'} and &tsmsg("[Msg] Setting reads information according to sam file [$samOrder][$fnSam]\n"); 
	$samOrder == 1 or $samOrder == 2 or &stopErr("[Err] samOrder must be 1 | 2 \n"); 
	my $arI = ($samOrder-1) * 2; 
	my %flag_aln_R1 = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,2=0,6=1,7=0;0=0,2=0' ) }; 
	my %flag_aln_R2 = %{ &SeqAlnSunhh::mk_flag( 'keep' => '0=1,2=0,6=0,7=1' ) }; 

	my $fhSam = &SeqAlnSunhh::openSam( $fnSam, $fmtSam, { 'wiH'=>0, 'exe_samtools' => $opts{'exe_samtools', 'verbose'=> $opts{'verbose'} } } ); 
	my %tmp_cnt = ('cntN_step'=>5e6); 
	while ( &wantLineC($fhSam) ) {
		&fileSunhh::log_section($., \%tmp_cnt) and &tsmsg("[Msg]   Processing $. line.\n"); 
		m/^\@/ and next; 
		my @ta = &splitL("\t", $_); 
		defined $rd_info_h->{'ID2idx'}{$ta[0]} or next; 
		if ( defined $flag_aln_R1{$ta[1]} ) {
			# This is a aligned R1 read. 
			my $cnt_mismat = &SeqAlnSunhh::cnt_sam_mismatch(\@ta, 'set_rna'); 
			my $newV = ($cnt_mismat == 0) ? 2 : 1; 
			&_ch_val_inAR( $rd_info_h->{ 'map_idx' }[ $rd_info_h->{'ID2idx'}{$ta[0]} ], $newV, $arI, 0 ); 
		} elsif ( defined $flag_aln_R2{$ta[1]} ) {
			my $cnt_mismat = &SeqAlnSunhh::cnt_sam_mismatch(\@ta, 'set_rna'); 
			my $newV = ($cnt_mismat == 0) ? 2 : 1; 
			&_ch_val_inAR( $rd_info_h->{ 'map_idx' }[ $rd_info_h->{'ID2idx'}{$ta[0]} ], $newV, $arI, 1 ); 
		} else {
			next; 
		}
	}
	return; 
}# sub set_info_bySam() 

sub _ch_val_inAR {
	my ( $ar, $v, $arI, $r12_add ) = @_; 
	$ar->[$arI+$r12_add] < $v and $ar->[$arI+$r12_add] = $v; 
	return; 
}# sub _ch_val_inAR() 

sub prepare_input {
	$opts{'exe_samtools'} //= 'samtools'; 
	$opts{'exe_perl'} //= 'perl'; 
	$opts{'pl_extractFq'} //= '/home/Sunhh/tools/github/NGS_data_processing/extract_fq_by_list.pl'; 
	$opts{'outPref'} //= 'out.'; 
	( defined $opts{'inSam1'} and defined $opts{'inSam2'} ) or &stopErr("[Err] Need -inSam1 and -inSam2\n"); 
	defined $opts{'inFq1'} or &stopErr("[Err] Need -inFq1\n"); 
	if (defined $opts{'fmtSam1'}) {
		$opts{'fmtSam1'} = lc($opts{'fmtSam1'}); 
		$opts{'fmtSam1'} =~ m/^(sam|bam)$/ or &stopErr("[Err] bad -fmtSam1 $opts{'fmtSam1'}\n"); 
	} else {
		$opts{'inSam1'} =~ m/\.(bam|sam)$/i and $opts{'fmtSam1'} = lc($1); 
	}
	if (defined $opts{'fmtSam2'}) {
		$opts{'fmtSam2'} = lc($opts{'fmtSam2'}); 
		$opts{'fmtSam2'} =~ m/^(sam|bam)$/ or &stopErr("[Err] bad -fmtSam2 $opts{'fmtSam2'}\n"); 
	} else {
		$opts{'inSam2'} =~ m/\.(bam|sam)$/i and $opts{'fmtSam2'} = lc($1); 
	}
	$opts{'fmtSam1'} //= 'sam'; 
	$opts{'fmtSam2'} //= 'sam'; 
	$opts{'verbose'} and &tsmsg("[Msg] Set : -fmtSam1 $opts{'fmtSam1'}   -fmtSam2 $opts{'fmtSam2'}\n"); 
	
	return; 
}# sub prepare_input() 

sub load_fqID {
	my ($fn) = @_; 
	$opts{'verbose'} and &tsmsg("[Msg] Loading fq file [$fn] start.\n"); 
	my @back_id; 
	my %back_hash; 
	my $fh = &openFH( $fn, '<' ) or &stopErr("[Err] Failed to open file [$fn]\n"); 
	my $cnt = -1; 
	while (<$fh>) {
		m/^\@(\S+)/ or &stopErr("[Err] Bad ID!\n"); 
		$cnt ++; 
		push(@back_id, [0, 0, 0, 0, $1]); # ( [rd1_idx1, rd2_idx1, rd1_idx2, rd2_idx2, rdID], [rd1_idx, rd2_idx], ... ) # _idx: 0 - unmapped , 1 - mapped with mismatch , 2 - mapped without mismatch 
		$back_hash{$1} = $cnt; 
		<$fh>; <$fh>; <$fh>; 
	}
	close($fh); 
	$opts{'verbose'} and &tsmsg("[Msg] Loading fq file [$fn] finished.\n"); 
	return({ 'ID2idx'=>\%back_hash, 'map_idx'=>[@back_id] }); 
}# load_fqID()

sub tell_info {

my $ii = <<II; 
# 0 : unmapped ; 1 : mapped with mismatch ; 2 : mapped with 0-mismatch 
# Separate input reads into XX groups : 
#   Group 1 : map with 0-mismatch in both SAMs; 
#   Group 2 : map with 0-mismatch in SAM1 and with mismatch in SAM2 ; 
#   Group 3 : map with 0-mismatch in SAM2 and with mismatch in SAM1 ; 
#   Group 4 : map with 0-mismatch in SAM1 and not mapped in SAM2 ; 
#   Group 5 : map with 0-mismatch in SAM2 and not mapped in SAM1 ; 
#   Group 6 : map with mismatch in both SAMs ; 
#   Group 7 : map with mismatch in SAM1 and not mapped in SAM2 ; 
#   Group 8 : map with mismatch in SAM2 and not mapped in SAM1 ; 
#   Group 9 : Others including - 'not mapped in any SAMs' ; 
# If input reads are paired, two ends are considered together. 
#   Mapping : both ends can be mapped ; 
#   0-mismatch : both ends are mapped with 0-mismatch ; 
# For mismatch calculation : 
#   Cigars [HSIDX] are all considered as mismatches. 
#   Because I use hisat2, so 'XM:i:n' is also considered. 
#
II
	&LogInforSunhh::usage($ii); 
}

