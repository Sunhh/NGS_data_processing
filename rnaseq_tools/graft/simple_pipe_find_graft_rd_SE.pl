#!/usr/bin/perl
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
	"src_fa:s@", # Source genome where transmitted reads come from; 
	"tgt_fa:s@", # Target/Current genome where transmitted reads don't belong to; 
	"in_fq1:s@",  # fastq_R1; 
		"max_mismat_cnt:i", # -1, Filter out alignments with high mismatch bases. -1 means no filter applied. 
		"max_mismat_ratio:f", # -1, Filter out alignments with high mismatch ratio [0-1). -1 means no filter aplied.
	"wrk_dir:s",        # './', Assign a directory to store all resulting bam files, instead of in current folder. 
	
	"doStep:s", 
	"pl_step1:s",  # Align rnaseq of target genome to src_fa and tgt_fa; 
	"pl_step2:s",  # Retrieve alignments of possible transmitted reads.
	"pl_step3:s",  # Check if possible transmitted reads are contained by non/self-grafted samples.
		"transRdBam:s@", # 
	"pl_stepX1:s", # find_transmit_stepX1.pl
	"pl_stepX2:s", # find_transmit_stepX2.pl

	"exe_samtools:s", 
	"exe_hisat2:s", 
	"pl_runHisat2_with2pass:s", 
	"pl_fix_NHnum:s", 
); 

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) }; 
my %gg; 
&setGlob(); 
&applyOpt(); 

$gg{'doStep'} eq 'all' and $gg{'doStep'} = '1,2'; 
for my $stepID (split/,/, $gg{'doStep'}) {
	$stepID =~ s!\s!!g; 
	if ( $stepID eq '1' ) {
		for my $tk (qw/in_fq1 src_fa tgt_fa/) {
			unless ( defined $opts{$tk}) {
				&tsmsg("[Err] -$tk is required.\n"); 
				&LogInforSunhh::usage($gg{'help_txt'}); 
			}
		}
	} elsif ( $stepID eq '2' ) {
		for my $tk (qw/src_fa tgt_fa/) {
			unless ( defined $opts{$tk}) {
				&tsmsg("[Err] -$tk is required.\n"); 
				&LogInforSunhh::usage($gg{'help_txt'}); 
			}
		}
	} elsif ( $stepID eq '3' ) {
		for my $tk (qw/bgRdBam/) {
			unless ( defined $opts{$tk}) {
				&tsmsg("[Err] -$tk is required.\n"); 
				&LogInforSunhh::usage($gg{'help_txt'}); 
			}
		}
	} elsif ( $stepID eq 'X1' ) {
		for my $tk (qw/src_fa in_fq1/) {
			unless ( defined $opts{$tk}) {
				&tsmsg("[Err] -$tk is required.\n"); 
				&LogInforSunhh::usage($gg{'help_txt'}); 
			}
		}
	} elsif ( $stepID eq 'X2' ) {
		; 
	}
}
for my $stepID (split(/,/, $gg{'doStep'})) {
	$stepID =~ s!\s!!g; 
	my $cmd = ''; 
	if ( $stepID eq '1' ) {
		# Step 1 : Align the reads to all references and combine the bam files to one. 
		# Input  : -pref , -in_fq1 , -src_fa , -tgt_fa and other parameters. 
		# Output : $gg{'pref'}_comb.bam
		-e $gg{'pl_step1'} or &stopErr("[Err] pl_step1 [$gg{'pl_step1'}] not exist.\n"); 
		$cmd = "perl $gg{'pl_step1'} "; 
		$cmd .= " -pref $gg{'pref'} "; 
		$cmd .= " -in_fq1 $gg{'in_fq1'}[0] "; 
		$cmd .= " -wrk_dir $gg{'wrk_dir'} "; 
		$cmd .= join(" ", " ", map { "-src_fa $_" } @{$gg{'src_fa'}}); 
		$cmd .= join(" ", " ", map { "-tgt_fa $_" } @{$gg{'tgt_fa'}}); 
		$cmd .= " -max_mismat_cnt $gg{'max_mismat_cnt'} "; 
		$cmd .= " -max_mismat_ratio $gg{'max_mismat_ratio'} "; 
		$cmd .= " -pl_runHisat2_with2pass $gg{'pl_runHisat2_with2pass'} "; 
		$cmd .= " -pl_fix_NHnum $gg{'pl_fix_NHnum'} "; 
		$cmd .= " -exe_samtools $gg{'exe_samtools'} "; 
		$cmd .= " -exe_hisat2   $gg{'exe_hisat2'} "; 
		&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	} elsif ( $stepID eq '2' ) {
		# Step 2 : Read in src_fa, tgt_fa and pref_comb.bam to count possible transmitted reads from source genome. 
		# Input  : -src_fa , -tgt_fa , -pref/-inBam 
		# Output : pref.src2tgt_rdStat ; pref.src2tgt_rdList
		-e $gg{'pl_step2'} or &stopErr("[Err] pl_step2 [$gg{'pl_step2'}] not exist.\n"); 
		$cmd = "perl $gg{'pl_step2'} "; 
		$cmd .= " -pref $gg{'pref'} "; 
		$cmd .= " -wrk_dir $gg{'wrk_dir'} "; 
		# $cmd .= " -inBam $gg{'wrk_dir'}/$gg{'pref'}_comb.bam "; 
		$cmd .= join(" ", " ", map { "-src_fa $_" } @{$gg{'src_fa'}}); 
		$cmd .= join(" ", " ", map { "-tgt_fa $_" } @{$gg{'tgt_fa'}}); 
		&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	} elsif ( $stepID eq '3' ) {
		# Step 3 : Check if possible transmitted reads are contained by non/self-grafted samples.
		# Input  : -wrk_dir , -pref ; -bgRdBam wrk_dir/pref.src2tgt_rd.bam 
		# Output : -outBam wrk_dir/pref.src2tgt_cleanRd.bam
		#           Here I add .step3 to avoid repeat. 
		$cmd = "perl $gg{'pl_step3'} "; 
		$cmd .= " -pref $gg{'pref'} -wrk_dir $gg{'wrk_dir'} "; 
		$cmd .= join(" ", " ", map { "-bgRdBam $_" } @{$gg{'bgRdBam'}}); 
		&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	} elsif ( $stepID eq 'X1' ) {
		# Step X1 : Prepare alignments from non/self-grafted samples to source genome; 
		#  Output : pref_comb.nonTrans.bam
		$cmd = "perl $gg{'pl_stepX1'} "; 
		$cmd .= " -pref $gg{'pref'} "; 
		$cmd .= join(" ", " ", map { "-src_fa $_" } @{$gg{'src_fa'}}); 
		$cmd .= join(" ", " ", map { "-in_fq1 $_" } @{$gg{'in_fq1'}}); 
		$cmd .= " -wrk_dir $gg{'wrk_dir'} "; 
		$cmd .= " -max_mismat_cnt $gg{'max_mismat_cnt'} "; 
		$cmd .= " -max_mismat_ratio $gg{'max_mismat_ratio'} "; 
		$cmd .= " -exe_samtools $gg{'exe_samtools'} "; 
		$cmd .= " -exe_hisat2   $gg{'exe_hisat2'} "; 
		$cmd .= " -pl_runHisat2_with2pass $gg{'pl_runHisat2_with2pass'} "; 
		&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	} elsif ( $stepID eq 'X2' ) {
		# Step X2 : Compress pref_comb.nonTrans.bam by removing duplicated reads. 
		#  Output : pref_comb.nonTransRmdup.bam
		$cmd = "perl $gg{'pl_stepX2'} "; 
		$cmd .= " -pref $gg{'pref'} "; 
		$cmd .= " -wrk_dir $gg{'wrk_dir'} "; 
		$cmd .= " -exe_samtools $gg{'exe_samtools'} "; 
		&exeCmd_1cmd($cmd) and &stopErr("[Err] Failed at cmd: $cmd\n"); 
	} else {
		&tsmsg("[Wrn] Skip unknown step [$stepID]\n"); 
	}
}


sub setGlob {
	$gg{'dir_abs'}      = &fileSunhh::_dirname( &fileSunhh::_abs_path($0) );
	$gg{'dir_link'}     = &fileSunhh::_dirname( &fileSunhh::_abs_path_4link($0) ); 
	$gg{'pref'}         = 'out'; 
	$gg{'in_fq1'}       = []; 
	$gg{'in_fq2'}       = 'NA'; 

	$gg{'exe_samtools'} = 'samtools'; 
	$gg{'exe_hisat2'}   = 'hisat2'; 
	$gg{'pl_runHisat2_with2pass'} = '/home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/runHisat2_with2pass.pl'; 
	$gg{'pl_fix_NHnum'}           = '/home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/fix_NHnum.pl'; 
	$gg{'max_mismat_cnt'}   = -1; 
	$gg{'max_mismat_ratio'} = -1; 
	$gg{'wrk_dir'}          = './'; 
	$gg{'doStep'}           = 'all'; 
	$gg{'pl_step1'}         = (-e "$gg{'dir_abs'}/find_transmit_step1.pl") ? "$gg{'dir_abs'}/find_transmit_step1.pl" : "$gg{'dir_link'}/find_transmit_step1.pl" ; 
	$gg{'pl_step2'}         = (-e "$gg{'dir_abs'}/find_transmit_step2.pl") ? "$gg{'dir_abs'}/find_transmit_step2.pl" : "$gg{'dir_link'}/find_transmit_step2.pl" ; 
	$gg{'pl_step3'}         = (-e "$gg{'dir_abs'}/find_transmit_step3.pl") ? "$gg{'dir_abs'}/find_transmit_step3.pl" : "$gg{'dir_link'}/find_transmit_step3.pl" ; 
	$gg{'pl_stepX1'}        = (-e "$gg{'dir_abs'}/find_transmit_stepX1.pl") ? "$gg{'dir_abs'}/find_transmit_stepX1.pl" : "$gg{'dir_link'}/find_transmit_stepX1.pl" ; 
	$gg{'pl_stepX2'}        = (-e "$gg{'dir_abs'}/find_transmit_stepX2.pl") ? "$gg{'dir_abs'}/find_transmit_stepX2.pl" : "$gg{'dir_link'}/find_transmit_stepX2.pl" ; 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# perl $0   -pref $gg{'pref'}   -src_fa fa   -tgt_fa fa   -in_fq1 in_RS1.fq.gz [-in_fq1 in_RS2.fq.gz]
#
#   -pref       [$gg{'pref'}] Output prefix 
#   -wrk_dir    [$gg{'wrk_dir'}] Directory storing result files. 
#   -src_fa     [faFile] \@ Indexed fasta files for source taxa where transmitted reads come from. 
#   -tgt_fa     [faFile] \@ Indexed fasta files for target taxa where transmitted reads go to. 
#   -in_fq1     [faFile] Sample fastq to be processed. 
#   
#   -max_mismat_cnt    [$gg{'max_mismat_cnt'}]   Parameter for pl_runHisat2_with2pass. -1 means no limitation. 
#   -max_mismat_ratio  [$gg{'max_mismat_ratio'}] Parameter for pl_runHisat2_with2pass. -1 means no limitation. 
#
#   -doStep     [$gg{'doStep'}] Could be '1,2' / 'X1,X2' ; 
#     : 1    -  Align the reads to all references and combine the bam files to one. 
#                 Produce : pref_comb.bam
#     : 2    -  Read in src_fa, tgt_fa and pref_comb.bam to count possible transmitted reads from source genome.
#                 Produce : pref.src2tgt_rdStat
#                           pref.src2tgt_rdList
#                           pref.src2tgt_rd.bam
#     # 3    -  Filter possible transmitted reads by comparison with nonTransRmdup.bam; 
#     : X1   -  Prepare alignments from non/self-grafted samples to source genome; 
#                 Produce : pref_comb.nonTrans.bam
#     : X2   -  Compress pref_comb.nonTrans.bam by removing duplicated reads. 
#                 Produce : pref_comb.nonTransRmdup.bam
#
#
#   -exe_samtools      [$gg{'exe_samtools'}] 
#   -exe_hisat2        [$gg{'exe_hisat2'}]
#   -pl_runHisat2_with2pass    [$gg{'pl_runHisat2_with2pass'}]
#   -pl_fix_NHnum              [$gg{'pl_fix_NHnum'}]
#
#   -pl_step1                  [$gg{'pl_step1'}]
#   -pl_step2                  [$gg{'pl_step2'}]
#   -pl_step3                  [$gg{'pl_step3'}]
#   -pl_stepX1                 [$gg{'pl_stepX1'}]
#   -pl_stepX2                 [$gg{'pl_stepX2'}]
################################################################################
HH

	return; 
}# setGlob() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	
	# For outer tools
	for my $k (qw/exe_samtools exe_hisat2 pl_runHisat2_with2pass pl_fix_NHnum/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	for my $k (qw/pref in_fq1 max_mismat_cnt max_mismat_ratio wrk_dir src_fa tgt_fa doStep pl_step1 pl_step2 pl_step3 pl_stepX1 pl_stepX2/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	return; 
}# applyOpt() 

