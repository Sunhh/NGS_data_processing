#!/usr/bin/perl
# 2016-08-11 The order of groups in comparison list matters, and now I want to use 2nd group as object, and 1st group as reference. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use wm97Sunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_snpTbl:s", # Acc131_mask.snp_addGmP 
	  "lis_chrID2num:s", # chrID \\t number 
	"in_wind:s",   # WM97_w10ks10k.wind 
	"in_annot:s",  # WM97_v6.annot.definition.annot
	"set_para:s@",  # 'med01_qtCut=0.03;med01_grpLen=5;med01_grpGood=3;med01_grpExt=0;'
	"firstAsObjPop!", 
	"chk_scripts!", 
	"use_sepRunXPCLR!", 
); 

my %glob; 
$glob{'set_para_default'}  = ['']; 
push(@{$glob{'set_para_default'}}, 'windTag       =w10ks10k; wind_size=10000; wind_step=10000;'); 
push(@{$glob{'set_para_default'}}, 'xpclr_w       =-w1 0.0005 100 100; xpclr_p=-p0 0.7;'); 
push(@{$glob{'set_para_default'}}, 'slct_sweep_wind.pl_slct_colN=5; slct_sweep_wind.pl_bpCnt_colN=4;'); 
push(@{$glob{'set_para_default'}}, 'med01_qtCut   =0.03;   med01_grpLen   =5;med01_grpGood=3;   med01_grpExt=0;'); 
push(@{$glob{'set_para_default'}}, 'med02_qtCut_01=0.20;med02_grpLen_01=2;med02_grpGood_01=1;med02_grpExt_01=0;'); 
push(@{$glob{'set_para_default'}}, 'med02_qtCut_02=0.10;med02_grpLen_02=1;med02_grpGood_02=1;med02_grpExt_02=0;'); 
push(@{$glob{'set_para_default'}}, 'med02_minWindLen=30000'); 
$opts{'set_para'} //= $glob{'set_para_default'}; 

$opts{'set_para_help'} = join('', map { "# $_\n" } @{$glob{'set_para_default'}}); 

for my $a1 ( @{$opts{'set_para'}} ) {
	$a1 =~ m!^[\s;]*$! and next; 
	for my $p1 (split(/;+/, $a1)) {
		$p1 =~ s!^\s+|\s+$!!g; 
		$p1 eq '' and next; 
		$p1 =~ m!^(\S+)\s*=\s*(\S.*?)$! or &stopErr("[Err] Failed to parse para [$p1]\n"); 
		my ($k, $v) = ($1, $2); 
		$glob{$k} = $v; 
	}
}
for my $a1 (@{$glob{'set_para_default'}}) {
	$a1 =~ m!^[\s;]*$! and next; 
	for my $p1 (split(/;+/, $a1)) {
		$p1 =~ s!^\s+|\s+$!!g; 
		$p1 eq '' and next; 
		$p1 =~ m!^(\S+)\s*=\s*(\S.*?)$! or &stopErr("[Err] Failed to parse para [$p1]\n"); 
		my ($k, $v) = ($1, $2); 
		$glob{$k} //= $v; 
	}
}


my $help_txt = <<HH; 
####################################################################################################
# perl $0 set01_Grp14_to_Grp13 out_dir
#
#   'set01_Grp14_to_Grp13' format : refPop must be in the front, and objPop with selection should be after refPop
#     sample_ID_1 \\t grpID_refPop
#     sample_ID_2 \\t grpID_refPop
#     sample_ID_3 \\t grpID_objPop
#     sample_ID_4 \\t grpID_objPop
#
# -help
#
# -firstAsObjPop   [Boolean] If given,   the 1st pop will be object population. 
#                            By default, the 2nd pop will be used as object pop. 
#
# -chk_scripts     [Boolean]  Only check if the scripts are available. 
#
# -in_snpTbl       [filename] Similar to file 'Acc131_mask.snp_addGmP', 
#                    Format : chr        \\t pos \\t cM                   \\t CLVAs_97103_GS1 \\t CLVAs_RZ-901_GS3
#                             WM97_Chr01 \\t 13  \\t 0.000104984264258213 \\t G               \\t G
# -in_wind         [filename] Similar to file 'WM97_w10ks10k.wind', 
#                    Format : chrID      \\t WindS \\t WindE \\t WindL \\t BpCnt
#                             WM97_Chr01 \\t 1     \\t 10000 \\t 10000 \\t 9914
#                             WM97_Chr01 \\t 10001 \\t 20000 \\t 10000 \\t 9953
# -in_annot        [filename] Similar to file 'WM97_v6.annot.definition.annot'
#                    Format : chrID       \\t start  \\t end    \\t strand \\t ....
#                             WM97_Chr00  \\t 258813 \\t 259526 \\t +      \\t Cla000003 \\t C14402000:7:720:+ \\t Mitochondrial transcription termination factor
#
# -lis_chrID2num   [filename] 
#                    Format : chrID       \\t chrNumber
#                             PG1         \\t 1
#                             scf_01      \\t 2
# -set_para        [string] Default is as following: 
# ==========================================================
$opts{'set_para_help'}
# ==========================================================
# 
#
#
####################################################################################################
HH

$opts{'chk_scripts'} and do { &set_scripts(); exit(1); }; 
!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

### Look for perl scripts; 
&set_scripts(); 

my $fn_list = shift; 
my $wrk_dir = shift; 
-e "$wrk_dir" and &stopErr("[Err] Existed $wrk_dir\n"); 
&tsmsg("[Rec] Begin [$0]\n"); 
mkdir($wrk_dir) or &stopErr("[Err] Failed to create directory [$wrk_dir]\n"); 

my $fn_snp_wiGmP = $opts{'in_snpTbl'}; 
my $fn_wind      = $opts{'in_wind'}; 
my $fn_annot     = $opts{'in_annot'}; 
my $fn_chrID2num = $opts{'lis_chrID2num'}; 
$fn_chrID2num //= ''; 
if ( defined $fn_chrID2num and $fn_chrID2num ne '' ) {
	my $fh = &openFH( $fn_chrID2num, '<' ); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		( defined $ta[0] and defined $ta[1] ) or &tsmsg("[Wrn] Skip chrID2num line: $_\n"); 
		defined $glob{'chr_id2num'}{$ta[0]} and &stopErr("[Err] Repeat chrID [$ta[0]]\n"); 
		defined $glob{'chr_num2id'}{$ta[1]} and &stopErr("[Err] Repeat number [$ta[1]]\n"); 
		$glob{'chr_id2num'}{$ta[0]} = $ta[1]; 
		$glob{'chr_num2id'}{$ta[1]} = $ta[0]; 
	}
	close ($fh); 
}

my %lis = %{ &load_comp_list( $fn_list ) }; 
# Here lis_A relates to genofile1 in XPCLR, which is used as object population. 
if ( $opts{'firstAsObjPop'} ) {
	&fileSunhh::write2file( "$wrk_dir/lis_A", join("\n", map { "$_\t$lis{'grpIDs'}[0]"; } @{$lis{'IDs'}{ $lis{'grpIDs'}[0] }})."\n" ); 
	&fileSunhh::write2file( "$wrk_dir/lis_B", join("\n", map { "$_\t$lis{'grpIDs'}[1]"; } @{$lis{'IDs'}{ $lis{'grpIDs'}[1] }})."\n" ); 
} else {
	&fileSunhh::write2file( "$wrk_dir/lis_A", join("\n", map { "$_\t$lis{'grpIDs'}[1]"; } @{$lis{'IDs'}{ $lis{'grpIDs'}[1] }})."\n" ); 
	&fileSunhh::write2file( "$wrk_dir/lis_B", join("\n", map { "$_\t$lis{'grpIDs'}[0]"; } @{$lis{'IDs'}{ $lis{'grpIDs'}[0] }})."\n" ); 
}

mkdir("$wrk_dir/input/"); 
&exeCmd_1cmd("$glob{'prepare_xpclr_input_wiGmP.pl'} $wrk_dir/input/input $fn_snp_wiGmP $wrk_dir/lis_A $wrk_dir/lis_B $fn_chrID2num") and &stopErr("[Err] here.\n"); 
$glob{'fh_o_wXPCLR'} = &openFH( "$wrk_dir/xpclr_$glob{'windTag'}" , '>' ); 
opendir DD,"$wrk_dir/input/" or &stopErr("[Err] Failed to opendir [$wrk_dir/input/]\n"); 
for my $d (readdir(DD)) {
	$d =~ m/^\./ and next; 
	$d =~ m!^input\.(\S+)\.snp! or next; 
	my $chrID = $1; 
	my $chrN; 
	if ( defined $glob{'chr_id2num'}{$chrID} ) {
		$chrN = $glob{'chr_id2num'}{$chrID}; 
	} else {
		$chrN = &wm97Sunhh::chrID_to_number( $chrID , 'WM97_Chr'); 
		$chrN =~ m!^\d+$! or &stopErr("[Err] Failed to convert chrID [$chrID] to number [$chrN]\n"); 
		defined $glob{'chr_num2id'}{$chrN} and &stopErr("[Err] Repeat chrN [$chrN] for differnt chrID [$glob{'chr_num2id'}{$chrN} $chrID]\n"); 
		$glob{'chr_id2num'}{$chrID} = $chrN; 
		$glob{'chr_num2id'}{$chrN}  = $chrID; 
	}
	my $i_pref = "$wrk_dir/input/input"; 
	if ( $opts{'use_sepRunXPCLR'} ) {
		&exeCmd_1cmd("$glob{'sep_run_xpclr.pl'} ${i_pref}.${chrID} ' $glob{'xpclr_w'} $chrN $glob{'xpclr_p'} '   300000   ${i_pref}.${chrID}.snp.out.xpclr.txt | grep -v process") and &stopErr("[Err] Stop sep_run_xpclr.pl\n"); 
	} else {
		&exeCmd_1cmd("XPCLR -xpclr ${i_pref}.${chrID}_g1.geno ${i_pref}.${chrID}_g2.geno ${i_pref}.${chrID}.snp ${i_pref}.${chrID}.snp.out -w1 0.0005 100 100 $chrN -p0 0.7 | grep -v process") and &stopErr("[Err] Failed to run XPCLR\n"); 
	}
	&exeCmd_1cmd("$glob{'cluster_xpclrscore.pl'} ${i_pref}.${chrID}.snp.out.xpclr.txt -wind_size $glob{'wind_size'} -wind_step $glob{'wind_step'} -wind_start 1 > ${i_pref}.${chrID}.snp.out.xpclr.txt.$glob{'windTag'}") and &stopErr("[Err] Stop cluster_xpclrscore.pl\n"); 
	my $fh = &openFH("${i_pref}.${chrID}.snp.out.xpclr.txt.$glob{'windTag'}", '<'); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		my $tk = "$ta[0]\t$ta[1]\t$ta[2]"; 
		defined $glob{'uniq_wind'}{$tk} and next; 
		$glob{'uniq_wind'}{$tk} = 1; 
		$ta[0] =~ m!^chrID$!i or do { defined $glob{'chr_num2id'}{$ta[0]} or &stopErr("[Err] Unknown chrN [$ta[0]]\n"); $ta[0] = $glob{'chr_num2id'}{$ta[0]}; }; 
		print { $glob{'fh_o_wXPCLR'} } join("\t", @ta)."\n"; 
	}
	close($fh); 
}
closedir(DD); 
close( $glob{'fh_o_wXPCLR'} ); 
&exeCmd_1cmd("$glob{'ColLink.pl'} $fn_wind -f1 $wrk_dir/xpclr_$glob{'windTag'} -keyC1 0,1 -keyC2 0,1 -add -Col1 4 -fill 'NA' > $wrk_dir/wind.xpclr.compare") and &stopErr("[Err] Stop $glob{'ColLink.pl'}\n"); 

# Method 01 
&exeCmd_1cmd("$glob{'slct_sweep_wind.pl'} -inRatioFile $wrk_dir/wind.xpclr.compare -qtCutoff $glob{'med01_qtCut'} -grpLen $glob{'med01_grpLen'} -grpGood $glob{'med01_grpGood'} -grpExtend $glob{'med01_grpExt'} -slct_colN $glob{'slct_sweep_wind.pl_slct_colN'} -bpCnt_colN $glob{'slct_sweep_wind.pl_bpCnt_colN'} -joinNeighbor > $wrk_dir/wind.xpclr.compare.med01_slct01") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("$glob{'ret_annot_by_loc.pl'}   -inLocLis $wrk_dir/wind.xpclr.compare.med01_slct01 -inAnnotLis $fn_annot   > $wrk_dir/wind.xpclr.compare.med01_slct01.annot") and &stopErr("[Err] Stop\n"); 

# Method 02 : PMID - 26358652 
&exeCmd_1cmd("$glob{'slct_sweep_wind.pl'} -inRatioFile $wrk_dir/wind.xpclr.compare -qtCutoff $glob{'med02_qtCut_01'} -grpLen $glob{'med02_grpLen_01'} -grpGood $glob{'med02_grpGood_01'} -grpExtend $glob{'med02_grpExt_01'} -slct_colN $glob{'slct_sweep_wind.pl_slct_colN'} -bpCnt_colN $glob{'slct_sweep_wind.pl_bpCnt_colN'} -joinNeighbor > $wrk_dir/wind.xpclr.compare.med02_slct01") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("$glob{'slct_sweep_wind.pl'} -inRatioFile $wrk_dir/wind.xpclr.compare.med02_slct01 -qtCutoff $glob{'med02_qtCut_02'} -grpLen $glob{'med02_grpLen_02'} -grpGood $glob{'med02_grpGood_02'} -grpExtend $glob{'med02_grpExt_02'} -slct_colN $glob{'slct_sweep_wind.pl_slct_colN'} -bpCnt_colN $glob{'slct_sweep_wind.pl_bpCnt_colN'} > $wrk_dir/wind.xpclr.compare.med02_slct02") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("awk ' NR == 1 || ( \$4 >= $glob{'med02_minWindLen'} )' $wrk_dir/wind.xpclr.compare.med02_slct02 > $wrk_dir/wind.xpclr.compare.med02_slct03") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("$glob{'ret_annot_by_loc.pl'}   -inLocLis $wrk_dir/wind.xpclr.compare.med02_slct03 -inAnnotLis $fn_annot   > $wrk_dir/wind.xpclr.compare.med02_slct03.annot") and &stopErr("[Err] Stop\n"); 

&tsmsg("[Rec] All done [$0]\n"); 

sub load_comp_list {
	my $fn = shift;
	my %back; 
	my $fh = &openFH($fn, '<'); 
	my %has; 
	while (&wantLineC($fh)) {
		my @ta=&splitL("\t", $_); 
		$back{'order'}{$ta[1]} //= $.; 
		defined $has{$ta[0]} and die "repeated $ta[0]\n"; 
		$has{$ta[0]} = 1; 
		push(@{$back{'IDs'}{$ta[1]}} , $ta[0]); 
	}
	close($fh); 
	$back{'grpIDs'} = [ sort { $back{'order'}{$a} <=> $back{'order'}{$b} } keys %{$back{'IDs'}} ]; 
	return(\%back); 
}

sub set_scripts {
	$glob{'bin_dir'}       //= &fileSunhh::_dirname( &fileSunhh::_abs_path($0) ); 
	&tsmsg("[Msg] bin_dir=[$glob{'bin_dir'}]\n"); 
	for my $pl (qw/ColLink.pl slct_sweep_wind.pl prepare_xpclr_input_wiGmP.pl sep_run_xpclr.pl cluster_xpclrscore.pl ret_annot_by_loc.pl /) {
		defined $glob{$pl} and next; 
		if ( -e "$glob{'bin_dir'}/$pl" ) {
			$glob{$pl} = "perl $glob{'bin_dir'}/$pl"; 
		} elsif ( my $ta = File::Which::which($pl) ) {
			$glob{$pl} = $ta; 
		} elsif ( -e "$glob{'bin_dir'}/../slct_sweep/$pl" ) {
			$glob{$pl} = "perl $glob{'bin_dir'}/../slct_sweep/$pl"; 
		} else {
			&stopErr("[Err] Failed to find script [$pl]\n"); 
		}
		&tsmsg("[Msg] Setting script [$pl] as [$glob{$pl}]\n"); 
	}
	return; 
}# set_scripts () 


