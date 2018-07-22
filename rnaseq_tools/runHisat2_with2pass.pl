#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use SeqAlnSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"in_pref_list:s", # Required. 
	"para_hisat2:s", 
	"max_mismat_cnt:i", "max_mismat_ratio:f", 
	"db_hisat2:s", 
	"exe_hisat2:s",   # hisat2
	"exe_samtools:s", # samtools
		"para_samtools_srt:s", 
	"wrk_dir:s",      # Work dir; 
	"cpuN:i", 
	"help!", 
); 

my %gg; 
&setGlob();

# Process single/paired fastq files (RGs) one by one; 
for my $q ( @{$gg{'fq_infor'}} ) {
	my %qh = %$q; 
	$gg{'MAX_PROCESSES'} = &LogInforSunhh::change_procN( $gg{'pm'}, $gg{'nprocF'}, $gg{'MAX_PROCESSES'} ); 
	my $pid = $gg{'pm'}->start and next; 
	### 
	my $oriDir = &fileSunhh::_abs_path("./"); 
	chdir($gg{'wrk_dir'}); 
	my $tmpDir = &fileSunhh::_abs_path_4link(&fileSunhh::new_tmp_dir('create' => 1)); 

	my $rgid = $qh{'RG'}; 
	my $rgLine = "--rg-id $rgid"; 
	$rgLine .= " --rg SM:$qh{'SM'}"; 
	$rgLine .= " --rg LB:$qh{'LB'}"; 
	$rgLine .= " --rg PL:$qh{'PL'}"; 
	$qh{'CN'} ne '' and $rgLine .= " --rg $qh{'CN'}"; 
	$qh{'DS'} ne '' and $rgLine .= " --rg $qh{'DS'}"; 
	$qh{'DT'} ne '' and $rgLine .= " --rg $qh{'DT'}"; 
	$qh{'PI'} ne '' and $rgLine .= " --rg $qh{'PI'}"; 
	$qh{'PG'} ne '' and $rgLine .= " --rg $qh{'PG'}"; 
	$qh{'PM'} ne '' and $rgLine .= " --rg $qh{'PM'}"; 

	

	my $cmd = ''; 

	# First pass : Align reads to db_index and generate splicesites with '--novel-splicesite-outfile'; 
	$cmd = ''; 
	$cmd .= "$gg{'exe_hisat2'} --novel-splicesite-outfile $tmpDir/splicesites.txt $gg{'para_hisat2'} -x $gg{'db_hisat2'} $rgLine "; 
	if ( $qh{'fq2'} eq '' ) {
		$cmd .= " -U $qh{'fq1'} "; 
	} else {
		$cmd .= " -1 $qh{'fq1'} -2 $qh{'fq2'} "; 
	}
	$cmd .= " -S /dev/null "; 
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed to run CMD: $cmd\n"); 
	# Second pass : Align reads to db_index with splicesites assigned by '--novel-splicesite-infile'; 
	$cmd = ''; 
	$cmd .= "$gg{'exe_hisat2'} --novel-splicesite-infile  $tmpDir/splicesites.txt $gg{'para_hisat2'} -x $gg{'db_hisat2'} $rgLine "; 
	if ( $qh{'fq2'} eq '' ) {
		$cmd .= " -U $qh{'fq1'} "; 
	} else {
		$cmd .= " -1 $qh{'fq1'} -2 $qh{'fq2'} "; 
	}
	$cmd .= " -S $tmpDir/$qh{'pref'}_x2.sam"; 
	&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed to run CMD: $cmd\n"); 

	# Filter bad alignments with high mismatch%; 
	# Convert sam file to bam file, and sort it. 
	&clean_sam2bam( "$tmpDir/$qh{'pref'}_x2.sam", "$qh{'pref'}_good.bam" ); 
	&exeCmd_1cmd("$gg{'exe_samtools'} sort $gg{'para_samtools_srt'} -o $qh{'pref'}_srt.bam $qh{'pref'}_good.bam") and &stopErr("[Err] Failed to sort bam file.\n"); 
	&fileSunhh::_rmtree("$qh{'pref'}_good.bam"); 

	&fileSunhh::_move( "$tmpDir/splicesites.txt", "$qh{'pref'}_splicesites.txt" ); 
	&fileSunhh::_rmtree("$qh{'pref'}_x2.sam"); 
	&fileSunhh::_rmtree($tmpDir); 

	chdir($oriDir); 
	&tsmsg("[Msg] Finished to align $qh{'pref'} and generated $qh{'pref'}_good.bam\n"); 
	$gg{'pm'}->finish; 
}
$gg{'pm'}->wait_all_children; 

# Clean the sam according to mismatches, and convert it to bam file. 
sub clean_sam2bam {
	my ($inSam, $outBam) = @_; 
	my %log_cnt; 
	%log_cnt = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
	open F,'<',"$inSam" or &stopErr("[Err] Failed to open sam file [$inSam]\n"); 
	open O,'|-',"$gg{'exe_samtools'} view -o $outBam -" or &stopErr("[Err] Failed to creat output bam [$outBam]\n"); 
	# select(O); $| = 1; 
	while (<F>) {
		&fileSunhh::log_section($., \%log_cnt) and &tsmsg("[Msg]   Reading $. line from file $inSam\n"); 
		m!^\@! and do { print O $_; next; }; 
		chomp; 
		my @ta = split(/\t/, $_); 
		my $statH = &stat_aln(\@ta); 
		if ( $gg{'max_mismat_ratio'} >= 0 ) {
			$statH->{'whole_mismat'} > $statH->{'read_len'} * $gg{'max_mismat_ratio'} and next; 
		}
		if ( $gg{'max_mismat_cnt'} >= 0 ) {
			$statH->{'whole_mismat'} > $gg{'max_mismat_cnt'} and next; 
		}
		print O "$_\n"; 
	}
	close(O); 
	close(F); 
	return; 
}# clean_sam2bam () 

sub stat_aln {
	my ($ar) = @_; 
	my %backH; 
	$backH{'cigarH'} = &SeqAlnSunhh::parseCigar( $ar->[5] ); 
	$backH{'read_len'} = $backH{'cigarH'}{'RdLen'}; 
	$backH{'whole_mismat'} = 0; 
	for my $tk (qw/Slen Hlen/) {
		defined $backH{'cigarH'}{$tk} and $backH{'whole_mismat'} += $backH{'cigarH'}{$tk}; 
	}
	for my $tb (@{$ar}[ 11 .. $#$ar ]) {
		$tb =~ m!^NM:i:(\d+)$! or next; 
		$backH{'whole_mismat'} += $1; 
		last; 
	}
	return(\%backH); 
}# stat_aln() 

# Return for paired: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>fq2, 'PL'=>'illumina', 'PU'=>RG}, {}, ... ])
# Return for single: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>'' , 'PL'=>'illumina', 'PU'=>RG}, {}, ... ])

sub setGlob {
	$gg{'exe_hisat2'} = 'hisat2'; 
	$gg{'exe_samtools'} = 'samtools'; 
	$gg{'para_samtools_srt'} = ' -@ 4 -m 5G '; 
	$gg{'para_hisat2'} = ' -p 4 --dta-cufflinks -q --phred33 '; 
	$gg{'max_mismat_ratio'} = -1; 
	$gg{'max_mismat_cnt'} = -1; 
	$gg{'oriDir'}     = &fileSunhh::_abs_path_4link("./"); 
	$gg{'wrk_dir'}    = &fileSunhh::_abs_path_4link("./"); 
	$gg{'cpuN'}       = 1; 
	
	$gg{'help_txt'} = <<HH; 
################################################################################
# perl $0  -in_pref_list pref_file  -db_hisat2 hisat2_index 
#
#  -help 
#  -cpuN                  [$gg{'cpuN'}] Number of RGs to be processed in parallel. 
#
#  -in_pref_list          [filename] Required. 
#                          I recommend to move the RGs with more reads to higher rank. 
#                          Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2 [ <tab> PL <tab> PU <tab> Others]
#                                   SM                RG_ID                 RG_LB              ReadOutPref      filename     filename  [    PL       PU       text]
#                                   Repeat            Repeat                Repeat             Unique           Unique       Unique    [    Rep      Rep      Rep]
#                                   GS109             GS109_Time1           GS109_Time1        GS109_Time1      W1_R1.fq.gz  W1_R2.fq.gz    illumina H0V2RADXX.1 NA
#                          Example of Others :
#                                   DS=GS109_WDM_Time1;WDM_L1_I121.R[12].clean.fastq.gz;Len_100_100;;DT=2014-02-01;;CN=BFC2013288
#                          'NA' means no input assigned.
#
#  -db_hisat2             [hisat2_index_path] Here I recommand to use index built without annotation. 
#                           This pipeline will align reads with 2-pass method used in hisat paper. 
#
#  -wrk_dir               ['./'] Assign a directory to store all resulting bam files, instead of in current folder. 
#  -para_hisat2           ['$gg{'para_hisat2'}']
#  -max_mismat_cnt        [$gg{'max_mismat_cnt'}] Filter out alignments with high mismatch bases. -1 means no filter applied. 
#  -max_mismat_ratio      [$gg{'max_mismat_ratio'}] Filter out alignments with high mismatch ratio [0-1). -1 means no filter aplied. 
#
#  -exe_hisat2            ['$gg{'exe_hisat2'}']
#  -exe_samtools          ['$gg{'exe_samtools'}']
#
#    -para_samtools_srt   ['$gg{'para_samtools_srt'}']
#
################################################################################
HH
	defined $opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	for my $t0 (qw/in_pref_list db_hisat2/) {
		defined $opts{$t0} or &stopErr("[Err] -$t0 required.\n"); 
	}
	$gg{'fq_infor'} = &load_prefList( $opts{'in_pref_list'} ); # Read in fastq information for bam output.
	$gg{'db_hisat2'} = &fileSunhh::_abs_path_4link( $opts{'db_hisat2'} ); 

	# replace default executable tools; 
	for my $exeTool (qw/exe_hisat2 exe_samtools/) {
		defined $opts{$exeTool} and $gg{$exeTool} = $opts{$exeTool}; 
		$gg{$exeTool} = &fileSunhh::_which( $gg{$exeTool} ); 
		defined $gg{$exeTool} or &stopErr("[Err] Failed to find -$exeTool [$opts{$exeTool}]\n"); 
	}
	# replace default parameters. 
	for my $pp (qw/para_hisat2 max_mismat_ratio max_mismat_cnt wrk_dir cpuN para_samtools_srt/) {
		defined $opts{$pp} and $gg{$pp} = $opts{$pp}; 
	}
	$gg{'wrk_dir'} = &fileSunhh::_abs_path_4link( $gg{'wrk_dir'} ); 
	-d $gg{'wrk_dir'} or do { mkdir($gg{'wrk_dir'}) or &stopErr("[Err] Failed to creat wrk_dir $opts{'wrk_dir'}\n"); }; 

	# Set multi-threads;
	$gg{'MAX_PROCESSES'} = $gg{'cpuN'}; 
	$gg{'nprocF'} = &fileSunhh::_abs_path_4link("$gg{'wrk_dir'}/general_Nproc"); 
	&fileSunhh::write2file($gg{'nprocF'}, "$gg{'MAX_PROCESSES'}\n", '>'); 
	$gg{'pm'} = &LogInforSunhh::get_pm( $gg{'MAX_PROCESSES'} ); 
	return; 
}# setGlob() 


# I would ignore lines heading with '#'
# Return for paired: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>fq2, 'PL'=>'illumina', 'PU'=>RG}, {}, ... ])
# Return for single: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>'' , 'PL'=>'illumina', 'PU'=>RG}, {}, ... ])
sub load_prefList {
	my $fn = shift;
	my $fh = &openFH($fn, '<');
	# Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2 [ <tab> PL <tab> PU <tab> Others]
	#          SM                RG_ID                 RG_LB              ReadOutPref      filename     filename  [    PL       PU       text]
	#          Repeat            Repeat                Repeat             Unique           Unique       Unique    [    Rep      Rep      Rep]
	#          GS109             GS109_Time1           GS109_Time1        GS109_Time1      W1_R1.fq.gz  W1_R2.fq.gz    illumina H0V2RADXX.1 NA
	# Example of Others :
	#          DS=GS109_WDM_Time1;WDM_L1_I121.R[12].clean.fastq.gz;Len_100_100;;DT=2014-02-01;;CN=BFC2013288
	my @back; 
	my %good_otherK = (
		'CN' => 1, 
		'DS' => 1, 
		'DT' => 1, 
		'PI' => 1, 
		'PG' => 1, 
		'PM' => 1
	); 
	while (<$fh>) {
		m/^\s*($|#)/ and next;
		chomp;
		my @ta = &splitL("\t", $_); 
		my ($sm, $rg, $lb, $pref, $fq1, $fq2, $pl, $pu, $others) = @ta;
		(defined $pl and $pl ne '' and $pl !~ m!^NA$!i) or $pl = 'illumina'; 
		(defined $pu and $pu ne '' and $pu !~ m!^NA$!i) or $pu = $rg; 
		my %th; # read file hash; 
		$fq2 //= '';
		($fq1 eq '' or $fq1 =~ m/^NA$/i) and do { &tsmsg("[Wrn] Skip bad line without fq1 [$_]\n"); next; };
		$fq2 =~ m/^NA$/i and $fq2 = ''; 
		$fq1 = &fileSunhh::_abs_path($fq1); 
		$fq2 ne '' and $fq2 = &fileSunhh::_abs_path($fq2); 
		$th{'SM'} = $sm;
		$th{'RG'} = $rg;
		$th{'LB'} = $lb;
		$th{'pref'} = $pref;
		$th{'fq1'} = $fq1;
		$th{'fq2'} = $fq2;
		$th{'PL'}  = $pl;
		$th{'PU'}  = $pu;
		if (defined $others and $others ne '' and $others !~ m!^NA$!i) {
			for my $t1 (split(/;;/, $others)) {
				$t1 =~ m!^([^\s=]+)=(.+)$! or &stopErr("[Err] Failed to parse Others text [$t1]\n"); 
				my ($k, $v) = ($1, $2); 
				defined $good_otherK{$k} or do { &tsmsg("[Wrn] Skip unknown RG tag [$k] for line: $_\n"); next; }; 
				$th{$k} = $v; 
			}
		}
		for my $k (keys %good_otherK) {
			$th{$k} //= ''; 
		}
		push(@back, \%th);
	}
	close($fn);
	return(\@back);
}# load_prefList()


