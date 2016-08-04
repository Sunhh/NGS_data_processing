#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use wm97Sunhh; 

!@ARGV and die "perl $0 set01_Grp14_to_Grp13 out_dir\n"; 

my $fn_list = shift; 
my $wrk_dir = shift; 
-e "$wrk_dir" and &stopErr("[Err] Existed $wrk_dir\n"); 
&tsmsg("[Rec] Begin [$0]\n"); 
mkdir($wrk_dir) or &stopErr("[Err] Failed to create directory [$wrk_dir]\n"); 


my $fn_snp_wiGmP = 'Acc131_mask.snp_addGmP'; 
my $fn_wind  = 'WM97_w10ks10k.wind'; 
my $fn_annot = 'WM97_v6.annot.definition.annot'; 

my %glob; 

$glob{'med01_qtCut'}   = 0.03; 
$glob{'med01_grpLen'}  = 5; 
$glob{'med01_grpGood'} = 3; 
$glob{'med01_grpExt'}  = 0; 


my %lis = %{ &load_comp_list( $fn_list ) }; 
&fileSunhh::write2file( "$wrk_dir/lis_A", join("\n", map { "$_\t$lis{'grpIDs'}[0]"; } @{$lis{'IDs'}{ $lis{'grpIDs'}[0] }})."\n" ); 
&fileSunhh::write2file( "$wrk_dir/lis_B", join("\n", map { "$_\t$lis{'grpIDs'}[1]"; } @{$lis{'IDs'}{ $lis{'grpIDs'}[1] }})."\n" ); 

mkdir("$wrk_dir/input/"); 
&exeCmd_1cmd("perl prepare_xpclr_input_wiGmP.pl $wrk_dir/input/input $fn_snp_wiGmP $wrk_dir/lis_A $wrk_dir/lis_B") and &stopErr("[Err] here.\n"); 
$glob{'fh_o_wXPCLR'} = &openFH( "$wrk_dir/xpclr_w10ks10k" , '>' ); 
opendir DD,"$wrk_dir/input/" or die; 
for my $d (readdir(DD)) {
	$d =~ m/^\./ and next; 
	$d =~ m!^input\.(\S+)\.snp! or next; 
	my $chrID = $1; 
	my $chrN = &wm97Sunhh::chrID_to_number( $chrID , 'WM97_Chr'); 
	my $i_pref = "$wrk_dir/input/input"; 
	&exeCmd_1cmd("perl sep_run_xpclr.pl   ${i_pref}.${chrID} ' -w1 0.0005 100 100 $chrN -p0 0.7 '   300000   ${i_pref}.${chrID}.snp.out.xpclr.txt | grep -v process") and &stopErr("[Err] Stop sep_run_xpclr.pl\n"); 
	&exeCmd_1cmd("perl cluster_xpclrscore.pl ${i_pref}.${chrID}.snp.out.xpclr.txt -wind_size 10000 -wind_step 10000 -wind_start 1 > ${i_pref}.${chrID}.snp.out.xpclr.txt.w10ks10k") and &stopErr("[Err] Stop cluster_xpclrscore.pl\n"); 
	my $fh = &openFH("${i_pref}.${chrID}.snp.out.xpclr.txt.w10ks10k", '<'); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		my $tk = "$ta[0]\t$ta[1]\t$ta[2]"; 
		defined $glob{'uniq_wind'}{$tk} and next; 
		$glob{'uniq_wind'}{$tk} = 1; 
		$ta[0] =~ m!^chrID$!i or $ta[0] = &wm97Sunhh::number_to_chrID( $ta[0] ); 
		print { $glob{'fh_o_wXPCLR'} } join("\t", @ta)."\n"; 
	}
	close($fh); 
}
closedir(DD); 
close( $glob{'fh_o_wXPCLR'} ); 
&exeCmd_1cmd("ColLink.pl $fn_wind -f1 $wrk_dir/xpclr_w10ks10k -keyC1 0,1 -keyC2 0,1 -add -Col1 4 -fill 'NA' > $wrk_dir/wind.xpclr.compare") and &stopErr("[Err] Stop ColLink.pl\n"); 

# Method 01 
&exeCmd_1cmd("perl slct_sweep_wind.pl -inRatioFile $wrk_dir/wind.xpclr.compare -qtCutoff $glob{'med01_qtCut'} -grpLen $glob{'med01_grpLen'} -grpGood $glob{'med01_grpGood'} -grpExtend $glob{'med01_grpExt'} -slct_colN 5 -bpCnt_colN 4 -joinNeighbor > $wrk_dir/wind.xpclr.compare.med01_slct01") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("perl ret_annot_by_loc.pl   -inLocLis $wrk_dir/wind.xpclr.compare.med01_slct01 -inAnnotLis $fn_annot   > $wrk_dir/wind.xpclr.compare.med01_slct01.annot") and &stopErr("[Err] Stop\n"); 

# Method 02 : PMID - 26358652 
&exeCmd_1cmd("perl slct_sweep_wind.pl -inRatioFile $wrk_dir/wind.xpclr.compare -qtCutoff 0.20 -grpLen 2 -grpGood 1 -grpExtend 0 -slct_colN 5 -bpCnt_colN 4 -joinNeighbor > $wrk_dir/wind.xpclr.compare.med02_slct01") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("perl slct_sweep_wind.pl -inRatioFile $wrk_dir/wind.xpclr.compare.med02_slct01 -qtCutoff 0.10 -grpLen 1 -grpGood 1 -grpExtend 0 -slct_colN 5 -bpCnt_colN 4 > $wrk_dir/wind.xpclr.compare.med02_slct02") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("awk ' NR == 1 || ( \$4 >= 3000 )' $wrk_dir/wind.xpclr.compare.med02_slct02 > $wrk_dir/wind.xpclr.compare.med02_slct03") and &stopErr("[Err] Stop\n"); 
&exeCmd_1cmd("perl ret_annot_by_loc.pl   -inLocLis $wrk_dir/wind.xpclr.compare.med02_slct03 -inAnnotLis $fn_annot   > $wrk_dir/wind.xpclr.compare.med02_slct03.annot") and &stopErr("[Err] Stop\n"); 

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


