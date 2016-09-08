#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Parallel::ForkManager; 
use mathSunhh; 
use fileSunhh; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"o_pref:s",      # Required ; 
	"grp_list:s",    # Format : individual_ID \\t anything; 
	"snp_tbl:s",     # Format : (Has header); chr \\t pos \\t indv_1 \\t indv_2 ... ; 
	"ncpu:i",        # 20 
	"wind_length:i", # 10000000 
	"wind_step:i",   # 1000000
	"pl_snpTbl_sepByWind:s", # /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/fst/snpTbl_sepByWind.pl 
	"pl_dealTbl:s",          # /home/Sunhh/tools/github/NGS_data_processing/deal_table.pl
	"jar_haploview:s",       # /data/Sunhh/src/Evolution/haploview/Haploview4.2.jar
	"maxNmissR:f",   # 0.40 
	"minMAF:f",      # 0.05 
	"tmp_dir:s",     # Default not assigned. 
	"chr_kln:s",     # in_chr.fa.key_len ; Format : chrID \\t chrLen \\t unique_number \\n 
	"maxdistance:i", # 1000 
	"hwcutoff:f",    # 0.001; 
	"memory:i",      # 10480 == 10G 
	"keep_tmp!",     
); 

my %glob;
my %default_para; 
$default_para{'ncpu'}        = 20; 
$default_para{'wind_length'} = 10e6; 
$default_para{'wind_step'}   = 1e6; 
$default_para{'pl_snpTbl_sepByWind'} = '/home/Sunhh/tools/github/NGS_data_processing/reseq_tools/fst/snpTbl_sepByWind.pl'; 
$default_para{'pl_dealTbl'}          = '/home/Sunhh/tools/github/NGS_data_processing/deal_table.pl'; 
$default_para{'jar_haploview'}       = '/data/Sunhh/src/Evolution/haploview/Haploview4.2.jar'; 
$default_para{'maxNmissR'}   = 0.10; 
$default_para{'minMAF'}      = 0.05; 
$default_para{'maxdistance'} = 1000; 
$default_para{'hwcutoff'}    = 0.001; 
$default_para{'memory'}      = 10480; # in Mb 
$default_para{'tmp_dir'}     = ''; 
 
&chk_para(); 
&set_para(); 
&run_pipe(); 
defined $glob{'keep_tmp'} or &fileSunhh::_rmtree( $glob{'tmp_dir'} ); 
&tsmsg("[Rec] All done. [$0]\n"); 

##### Sub-routines. 
sub set_para {
	for my $tk (keys %opts) {
		$glob{$tk} = $opts{$tk}; 
	}

	for my $tk (keys %default_para) {
		$glob{$tk} //= $default_para{$tk}; 
	}

	if ( defined $glob{'tmp_dir'} and $glob{'tmp_dir'} ne '' ) {
		-d $glob{'tmp_dir'} and &stopErr("[Err] tmp_dir [$glob{'tmp_dir'}] already exist.\n"); 
	} else {
		$glob{'tmp_dir'} = fileSunhh::new_tmp_dir(); 
		defined $glob{'tmp_dir'} or &stopErr("[Err] Failed to find temporary directory.\n"); 
	}
	mkdir($glob{'tmp_dir'}); 

	$glob{'pm'} = new Parallel::ForkManager( $glob{'ncpu'} ); 

	return; 
}#set_para() 

sub chk_para {
	my $help_txt = <<HH; 
####################################################################################################
perl $0 -snp_tbl in_snp.tbl -grp_list indiv_grp.list -o_pref out_prefix -chr_kln in_chr.fa.key_len_number

-help 

-snp_tbl           [] Format: (Has header); chr \\t pos \\t indv_1 \\t indv_2 ... ; 
-grp_list          [] Format: individual_ID \\t anything; 
-chr_kln           [] Format: chrID \\t chrLen \\t unique_number \\n


-ncpu              [$default_para{'ncpu'}]
-wind_length       [$default_para{'wind_length'}]
-wind_step         [$default_para{'wind_step'}]
-maxNmissR         [0-1] Maximum missing ratio allowed. Default is [$default_para{'maxNmissR'}]
-minMAF            [0-1] Minimum missing ratio allowed. Default is [$default_para{'minMAF'}]

For haploview parameters : 
-maxdistance       [$default_para{'maxdistance'}]
-hwcutoff          [$default_para{'hwcutoff'}]
-memory            [$default_para{'memory'}]


-pl_snpTbl_sepByWind      [$default_para{'pl_snpTbl_sepByWind'}]
-pl_dealTbl               [$default_para{'pl_dealTbl'}]
-jar_haploview            [$default_para{'jar_haploview'}]

-tmp_dir                  [$default_para{'tmp_dir'}]
-keep_tmp                 [Boolean] If given, '-tmp_dir' won't be removed. 

...

####################################################################################################
HH

	defined $opts{'help'} and &LogInforSunhh::usage($help_txt); 
	for (qw/snp_tbl o_pref grp_list chr_kln/) {
		defined $opts{$_} or &LogInforSunhh::usage($help_txt); 
	}
	return; 
}# chk_para() 

sub filter_site {
	my ($in_fn, $out_fn, $max_missR, $min_maf) = @_; 
	$max_missR //= $glob{'maxNmissR'}; 
	$min_maf   //= $glob{'minMAF'}; 

	my %dna2arr; 
	my %dna2d; 
	for my $g1 (qw/A T G C a t g c/) {
		my $g1_uc = uc($g1); 
		for my $g2 (qw/A T G C a t g c/) {
			my $g2_uc = uc($g2); 
			my @arr = sort ( $g1_uc, $g2_uc ); 
			$dna2arr{"$g1$g2"} //= [@arr]; 
			my $hete_uc = uc( SNP_tbl::dna_b2d("$g1_uc$g2_uc") ); 
			my $hete_lc = lc( $hete_uc ); 
			$dna2arr{$hete_uc} //= [@arr]; 
			$dna2arr{$hete_lc} //= [@arr]; 
			$dna2d{"$g1$g2"}   //= $hete_uc; 
			$dna2d{$hete_uc}   //= $hete_uc; 
		}
		$dna2arr{$g1} //= [ $g1_uc , $g1_uc ]; 
		$dna2d{$g1}   //= $g1_uc; 
	}

	my $in_fh  = &openFH($in_fn , '<'); 
	my $out_fh = &openFH($out_fn, '>'); 
	my $has_siteN = 0; 
	SNP_LINE: 
	while (<$in_fh>) {
		$. == 1 and do { print {$out_fh} $_; next; }; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		my %cnt; 
		$cnt{'allN'} = $#ta-2+1; 
		for my $tb (@ta[2 .. $#ta]) {
			$tb = uc($tb); 
			$tb eq 'N' and do { $cnt{'missN'} ++; next; }; 
			defined $dna2arr{$tb} or next SNP_LINE; 
			defined $dna2d{$tb} or &stopErr("[Err] undefined genotype [$tb]\n"); 
			$tb = $dna2d{$tb}; 
			$cnt{'alleleN'}{ $dna2arr{$tb}[0] } ++; 
			$cnt{'alleleN'}{ $dna2arr{$tb}[1] } ++; 
			$cnt{'alleleAN'} += 2; 
		}
		$cnt{'alleleTN'} = [ map { [$_ , $cnt{'alleleN'}{$_}] } keys %{$cnt{'alleleN'}} ]; 
		@{ $cnt{'alleleTN'} } == 2 or next SNP_LINE; 
		$cnt{'missN'} //= 0; 
		$cnt{'missN'} <= $cnt{'allN'} * $max_missR or next SNP_LINE; 
		$cnt{'alleleTN'}[0][1] >= $cnt{'alleleAN'} * $min_maf or next SNP_LINE; 
		$cnt{'alleleTN'}[1][1] >= $cnt{'alleleAN'} * $min_maf or next SNP_LINE; 
		print {$out_fh} join("\t", @ta)."\n"; 
		$has_siteN ++; 
	}
	close($out_fh); 
	close($in_fh); 
	return($has_siteN); 
}# filter_site() 

sub cnvt_snpTbl_to_ped {
	my ($in_fn, $ped_fn, $inf_fn) = @_; 
	
	my %dna2arr; 
	for my $g1 (qw/A T G C/) {
		my $g1_uc = uc($g1); 
		for my $g2 (qw/A T G C/) {
			my $g2_uc = uc($g2); 
			my @arr = sort ( $g1_uc, $g2_uc ); 
			$dna2arr{"$g1$g2"} //= [@arr]; 
			my $hete_uc = uc( SNP_tbl::dna_b2d("$g1_uc$g2_uc") ); 
			my $hete_lc = lc( $hete_uc ); 
			$dna2arr{$hete_uc} //= [@arr]; 
			$dna2arr{$hete_lc} //= [@arr]; 
		}
		$dna2arr{$g1} //= [ $g1_uc , $g1_uc ]; 
	}
	$dna2arr{'N'} = ['N', 'N']; 
	my %a2n = qw( A 1 C 2 G 3 T 4 N 0 ); 
	my %dna2arrTxt; 
	for my $tk (keys %dna2arr) {
		$dna2arrTxt{$tk} = $a2n{ $dna2arr{$tk}[0] } . " " . $a2n{ $dna2arr{$tk}[1] } ; 
	}

	my $in_fh  = &openFH( $in_fn , '<' ); 
	my $ped_fh = &openFH( $ped_fn, '>' ); 
	my $inf_fh = &openFH( $inf_fn, '>' ); 
	my @ped_line; 
	my $snp_idx = 0; 
	while (<$in_fh>) {
		$. == 1 and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		
		my %cnt; 
		for (my $i=2; $i<@ta; $i++) {
			my $tb = $ta[$i]; 
			defined $dna2arr{$tb} or &stopErr("[Err] Bad genotype [$tb]\n"); 
			$ped_line[$i] .= "\t$dna2arrTxt{$tb}"; 
			if ( $tb ne 'N' ) {
				$cnt{'al_N'}{ $dna2arr{$tb}[0] } ++; 
				$cnt{'al_N'}{ $dna2arr{$tb}[1] } ++; 
				$cnt{'al_NN'} += 2; 
			}
		}
		$cnt{'al_genoN'} = [ map { [ $_, $cnt{'al_N'}{$_} ] } sort { $cnt{'al_N'}{$b} <=> $cnt{'al_N'}{$a} || $a cmp $b } keys %{$cnt{'al_N'}} ]; 
		@{$cnt{'al_genoN'}} == 2 or &stopErr("[Err] Bad line for ped : $_\n"); 
		$snp_idx ++; 
		print {$inf_fh} join("\t", $snp_idx, $ta[1], $cnt{'al_genoN'}[0][0], sprintf("%.4f", $cnt{'al_genoN'}[0][1]/$cnt{al_NN}), $cnt{'al_genoN'}[1][0], sprintf("%.4f", $cnt{'al_genoN'}[1][1]/$cnt{al_NN}))."\n"; 
	}
	for ( my $i=2; $i<@ped_line; $i++ ) {
		print {$ped_fh} "$i\t$i\t0\t0\t0\t0" . $ped_line[$i] . "\n"; 
	}
	close($ped_fh); 
	close($inf_fh); 
	close($in_fh); 
}# cnvt_snpTbl_to_ped() 

sub bin_LD {
	my ($ld_fn, $prev_wSE) = @_; 
	$prev_wSE //= [ -1, -1, '.info' ]; 
	my $has_overlap = 1; 
	$prev_wSE->[0] > 0 or $has_overlap = 0; 
	my %idx2pos; 
	if ( $has_overlap == 1 ) {
		%idx2pos = map { $_->[0] => $_->[1] } &fileSunhh::load_tabFile( $prev_wSE->[2] ); 
	}

	my $ld_fh  = &openFH($ld_fn,'<'); 
	my %linked_2nd; 
	my %cnt; 
	while (<$ld_fh>) {
		$. == 1 and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		defined $linked_2nd{$ta[0]} and next; 
		if ($ta[4] >= 0.9) {
			$linked_2nd{$ta[1]} = 1; 
		}

		# Check if we have overlapping pairs with the previous one. 
		if ( $has_overlap == 1 ) {
			if ( $idx2pos{$ta[0]} <= $prev_wSE->[1] ) {
				$idx2pos{$ta[1]} <= $prev_wSE->[1] and next; 
			} else {
				$has_overlap = 0; 
			}
		}

		$cnt{$ta[7]}{'pairN'} ++; 
		$cnt{$ta[7]}{'sum_r2'} += $ta[4]; 
		$cnt{$ta[7]}{'sum_Dp'} += $ta[2]; 
	}
	my $bin_fh = &openFH("$ld_fn.cnt", '>'); 
	for my $tk (sort { $a <=> $b } keys %cnt) {
		print {$bin_fh} join("\t", $tk, $cnt{$tk}{'pairN'}, $cnt{$tk}{'sum_r2'}, $cnt{$tk}{'sum_Dp'}, $cnt{$tk}{'sum_r2'}/$cnt{$tk}{'pairN'}, $cnt{$tk}{'sum_Dp'}/$cnt{$tk}{'pairN'})."\n"; 
	}
	close($bin_fh); 
	close($ld_fh); 
	return; 
}# bin_LD() 

sub run_pipe {
	# 1. selecte individuals; 
	$glob{'k2ln'} = { map { $_->[0] => [$_->[1], $_->[2]] } &fileSunhh::load_tabFile( $glob{'chr_kln'} ) }; 
	&exeCmd_1cmd("perl $glob{'pl_dealTbl'} $glob{'snp_tbl'} -colByTbl $glob{'grp_list'} -colByTbl_also 0,1 > $glob{'tmp_dir'}/basic.snp") and &stopErr("[Err] Stop at selecting indv\n"); 
	# 2. split windows; in "$glob{'tmp_dir'}/sep_wind/" ; 
	&exeCmd_1cmd("perl $glob{'pl_snpTbl_sepByWind'} -tmp_dir $glob{'tmp_dir'}/sep_wind -snp_tbl $glob{'tmp_dir'}/basic.snp -skipSort -ncpu $glob{'ncpu'} -out $glob{'tmp_dir'}/wind_list -wind_start 1 -wind_end_useMax -wind_length $glob{'wind_length'} -wind_step $glob{'wind_step'} -chr_colN 0 -pos_colN 1 -skipHN 1") and &stopErr("[Err] Failed to separate windows.\n"); 
	my @wind_list = &fileSunhh::load_tabFile( "$glob{'tmp_dir'}/wind_list" , 1 ); # ( [ wind_1_fn, chrID, windS, windE, windLen ], [ ... ], ... ) 
	# 3. process each window; 
	for my $a1 ( @wind_list ) {
		my $pid = $glob{'pm'}->start and next; 

		# a1. filter missing and minMAF, and allele number by site; 
		my ($wind_fn, $wind_chrID, $wind_S, $wind_E) = @$a1; 
		defined $glob{'k2ln'}{$wind_chrID} or next; 
		my $wind_S_prev = $wind_S - $glob{'wind_step'}; 
		my $wind_E_prev = $wind_S_prev + $glob{'wind_length'} - 1; 
		my $has_siteN = &filter_site( $wind_fn, "$wind_fn.filter", $glob{'maxNmissR'}, $glob{'minMAF'} ); 
		if ( $has_siteN > 0 ) {
			# a2. generate .ped and .info files. 
			&cnvt_snpTbl_to_ped( "$wind_fn.filter", "$wind_fn.ped", "$wind_fn.info" ); 
			# a3. run haploview
			&exeCmd_1cmd("java -jar $glob{'jar_haploview'} -n -pedfile $wind_fn.ped -info $wind_fn.info -log $wind_fn.log -dprime -maxdistance $glob{'maxdistance'} -minMAF $glob{'minMAF'} -hwcutoff $glob{'hwcutoff'} -memory $glob{'memory'}") and &stopErr("[Err] Stop to run haploview for $wind_fn\n"); 
			# a4. summarize .LD file; 
			&bin_LD( "$wind_fn.ped.LD" , [$wind_S_prev, $wind_E_prev, "$wind_fn.info" ] ); # Write $wind_fn.ped.LD.cnt, format : v_Dist \\t pairs_of_dist \\t sum_R2 \\t sum_D' \\t avg_R2 \\t avg_D' ;
		}
		$glob{'pm'}->finish; 
	}
	$glob{'pm'}->wait_all_children; 
	# 4. collect .ped.LD.cnt files 
	my %cnt; 
	for my $a1 ( @wind_list ) {
		my ($wind_fn, $wind_chrID, $wind_S, $wind_E) = @$a1; 
		-e "$wind_fn.ped.LD.cnt" or do { &tsmsg("[Wrn] Missing file [$wind_fn.ped.LD.cnt]\n"); next; }; 
		my $c_fh = &openFH("$wind_fn.ped.LD.cnt", '<'); 
		while (<$c_fh>) {
			chomp; 
			my @ta = &splitL("\t", $_); 
			$cnt{$ta[0]}{'pairN'}  += $ta[1]; 
			$cnt{$ta[0]}{'sum_r2'} += $ta[2]; 
			$cnt{$ta[0]}{'sum_Dp'} += $ta[3]; 
		}
		close($c_fh); 
	}
	my $out_c_fh = &openFH( "$glob{'o_pref'}.LD_cnt" , '>'); 
	print {$out_c_fh} join("\t", qw/Dist pairNum sum_R2 sum_Dp Avg_r2 Avg_Dp/)."\n"; 
	for my $v_dist (sort {$a <=> $b} keys %cnt) {
		print {$out_c_fh} join("\t", 
		  $v_dist, 
		  $cnt{$v_dist}{'pairN'}, 
		  $cnt{$v_dist}{'sum_r2'}, 
		  $cnt{$v_dist}{'sum_Dp'}, 
		  $cnt{$v_dist}{'sum_r2'}/$cnt{$v_dist}{'pairN'}, 
		  $cnt{$v_dist}{'sum_Dp'}/$cnt{$v_dist}{'pairN'}
		)."\n"; 
	}
	close($out_c_fh); 
}# run_pipe() 


