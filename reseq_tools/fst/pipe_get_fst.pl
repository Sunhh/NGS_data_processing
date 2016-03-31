#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Parallel::ForkManager; 
use mathSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"ind2grp_list:s", "snp_tbl:s", 
	"o_pref:s", 
	"ncpu:i", # 20 
	"wind_length:i", "wind_step:i", 
	"pl_snpTbl_sepByWind:s", 
	"pl_cnvt_tbl2fstat:s", 
	"pl_run_hierfstat:s", "maxNmissR:f", "minGrp1N:i", "minGrp2N:i", "rmNegNeiFst!", "rmNegWcFst!", 
	"pl_join_fst_siteChrPos:s", 
	"exe_Rscript:s", 
); 

my %para; 

&set_opts(); 

&chk_opts(); 

&run_pipe(); 

&tsmsg("[Rec] All done for $0\n"); 

sub run_pipe {
	&exeCmd_1cmd("perl $opts{'pl_snpTbl_sepByWind'} -snp_tbl $opts{'snp_tbl'} -skipSort -ncpu $opts{'ncpu'} -out $opts{'o_pref'}.wind_list -wind_start 1 -wind_end_useMax -wind_length $opts{'wind_length'} -wind_step $opts{'wind_step'} -chr_colN 0 -pos_colN 1 -skipHN 1"); 
	my ($all_fstInListFn, $wid2cp_href) = &fmt_wlist( "$opts{'o_pref'}.wind_list" ); 
	my $dvd_fstInList = &sep_list( $all_fstInListFn, $opts{'ncpu'} ); 

	my $pm = new Parallel::ForkManager( $opts{'ncpu'} ); 
	for my $fn (@$dvd_fstInList) {
		my $pid = $pm->start and next; 
		open F,'<',"$fn" or die; 
		while (<F>) {
			chomp; 
			my @ta = split(/\t/, $_);
			&exeCmd_1cmd("perl $opts{'pl_cnvt_tbl2fstat'} $ta[0] -ind2grp_list $opts{'ind2grp_list'} -o_mrk_info $ta[2] > $ta[0].fmt"); 
			&exeCmd_1cmd("mv $ta[0].fmt $ta[0]"); 
		}
		close F; 
		&exeCmd_1cmd("perl $opts{'pl_run_hierfstat'} $para{'maxNmissR'} $para{'minGrp1N'} $para{'minGrp2N'} $para{'rmNegNeiFst'} $para{'rmNegWcFst'} -fst_in $fn -inList -exe_Rscript $opts{'exe_Rscript'} ") and &stopErr("[Err] Failed to execute $opts{'pl_run_hierfstat'}\n"); 
		$pm->finish; 
	}
	$pm->wait_all_children; 

	&exeCmd_1cmd("perl $opts{'pl_join_fst_siteChrPos'} $all_fstInListFn -addSuff .fst.perSiteChrPos > $opts{'o_pref'}.fst.perSiteChrPos"); 
	&exeCmd_1cmd("perl $opts{'pl_join_fst_siteChrPos'} $all_fstInListFn -perWind -addSuff .fst.perWindLine > $opts{'o_pref'}.fst.perWindLine"); 

	&addChrPos( "$opts{'o_pref'}.fst.perWindLine", $wid2cp_href, "$opts{'o_pref'}.fst.perWindChrPos" ); 

	&cleanup([@$dvd_fstInList, "$opts{'o_pref'}.fst.perWindLine", "$opts{'o_pref'}.wind_list", "$opts{'o_pref'}.wind_list.fmt.list"], []); # There is a sub_dir left, but I don't have directly data recording it. 

	return; 
}

sub cleanup {
	# $_[0] : [files]
	# $_[1] : [dirs]
	&tsmsg("[Msg] Clean temporary files.\n"); 
	for (@{$_[0]}) {
		unlink($_); 
	}
	defined $_[1] or return 0; 
	@{$_[1]} > 0 or return 0; 
	&fileSunhh::_rmtree($_[1]); 
	return 0; 
}

sub addChrPos {
	# $_[0]
	# $_[1] 
	# $_[2] 
	my $ifh = &openFH($_[0], '<'); 
	my $ofh = &openFH($_[2], '>'); 
	while (<$ifh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($. == 1) {
			print {$ofh} join("\t", qw/chr windS windE windL/, @ta)."\n"; 
			next; 
		}
		print {$ofh} join("\t", @{$_[1]->{$ta[0]}}[0 .. 3], @ta)."\n"; 
	}
	close($ofh); 
	close($ifh); 
	return; 
}


sub fmt_wlist {
	# $_[0] : $fn
	my %wid_to_chrPos; 
	my $fh = &openFH($_[0],'<'); 
	my $ofile = "$_[0].fmt.list"; 
	my $ofh = &openFH($ofile,'>'); 
	while (<$fh>) {
		chomp; 
		my @ta=split(/\t/, $_); 
		print {$ofh} join("\t", "$ta[0]", "$ta[0]", "$ta[0].mrkInf", @ta[1..$#ta])."\n"; 
		$wid_to_chrPos{$ta[0]} = [ @ta[ 1 .. $#ta ] ]; 
	}
	close($ofh); 
	close($fh); 
	return ($ofile, \%wid_to_chrPos); 
}

sub sep_list {
	# $_[0] : in.wind_list 
	# $_[1] : ncpu 
	my $fh = &openFH($_[0], '<'); 
	my @dd; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		push(@dd, $_); 
	}
	close($fh); 
	my $dvd_dd = &mathSunhh::dvd_array(\@dd, $_[1], 1); 
	my @back_ff; 
	for my $a1 (@$dvd_dd) {
		my $fn = &fileSunhh::new_tmp_file(); 
		my $ofh = &openFH($fn, '>'); 
		push(@back_ff, $fn); 
		for (@$a1) {
			print {$ofh} join('', $_)."\n"; 
		}
		close($ofh); 
	}
	return \@back_ff; 
}


sub chk_opts {
	my $help_txt = <<HH; 

perl $0 -snp_tbl in_snp.tbl -ind2grp_list indiv_to_grpNum -o_pref out_prefix

-help 

-ncpu              [$opts{'ncpu'}]
-wind_length       [$opts{'wind_length'}]
-wind_step         [$opts{'wind_step'}]
-maxNmissR         [0-1] Maximum missing ratio allowed. Default is no control. 
-rmNegNeiFst       [Bool] Set NA to sites with negative Nei_Fst if given. 
-rmNegWcFst        [Bool] Set NA to sites with negative Wc_Fst  if given. 
...

HH

	defined $opts{'help'} and &LogInforSunhh::usage($help_txt); 
	for (qw/snp_tbl o_pref ind2grp_list/) {
		defined $opts{$_} or &LogInforSunhh::usage($help_txt); 
	}
	return; 
}

sub set_opts {
	$opts{'pl_snpTbl_sepByWind'} //= '/home/Sunhh/tools/github/NGS_data_processing/reseq_tools/fst/snpTbl_sepByWind.pl'; 
	$opts{'pl_cnvt_tbl2fstat'} //= '/home/Sunhh/tools/github/NGS_data_processing/reseq_tools/fst/cnvt_tbl2fstat.pl'; 
	$opts{'pl_run_hierfstat'} //= '/home/Sunhh/tools/github/NGS_data_processing/reseq_tools/fst/run_hierfstat.pl'; 
	$opts{'pl_join_fst_siteChrPos'} //= '/home/Sunhh/tools/github/NGS_data_processing/reseq_tools/fst/join_fst_siteChrPos.pl'; 
	$opts{'exe_Rscript'} //= '~/bin/Rscript'; 

	$opts{'ncpu'} //= 20; 
	$opts{'wind_length'} //= 10000; 
	$opts{'wind_step'} //= $opts{'wind_length'}; 

	$para{'maxNmissR'}   = ( defined $opts{'maxNmissR'} ) ? " -maxNmissR $opts{'maxNmissR'} " : "" ; 
	$para{'minGrp1N'}    = ( defined $opts{'minGrp1N'} )  ? " -minGrp1N $opts{'minGrp1N'} "   : "" ; 
	$para{'minGrp2N'}    = ( defined $opts{'minGrp2N'} )  ? " -minGrp2N $opts{'minGrp2N'} "   : "" ; 
	$para{'rmNegNeiFst'} = ( $opts{'rmNegNeiFst'} ) ? " -rmNegNeiFst " : ""; 
	$para{'rmNegWcFst'}  = ( $opts{'rmNegWcFst'}  ) ? " -rmNegWcFst  " : ""; 

	return; 
}


