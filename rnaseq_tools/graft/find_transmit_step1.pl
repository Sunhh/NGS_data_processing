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
	"in_fq1:s",  # fastq_R1; 
	"in_fq2:s",  # fastq_R2; Not used yet. 
		"max_mismat_cnt:i", # -1, Filter out alignments with high mismatch bases. -1 means no filter applied. 
		"max_mismat_ratio:f", # -1, Filter out alignments with high mismatch ratio [0-1). -1 means no filter aplied.
	"wrk_dir:s",        # './', Assign a directory to store all resulting bam files, instead of in current folder. 
	
	"outBam:s",  # The combined bam file output. This will replace the output file "$gg{'wrk_dir'}/$gg{'pref'}_comb.bam"

	"exe_samtools:s", 
	"exe_hisat2:s", 
	"pl_runHisat2_with2pass:s", 
	"pl_fix_NHnum:s", 
); 

my %flag_UN = %{ &SeqAlnSunhh::mk_flag( 'keep' => '2=1' ) }; 
my %gg; 
&setGlob(); 
&applyOpt(); 
&step1_run_aln(); 

sub step1_run_aln {
	my $cmd = ''; 
	# &fileSunhh::write2file("$gg{'wrk_dir'}/$gg{'pref'}.pref_RS", '', '>'); 
	# Produce : 
	#   $gg{'wrk_dir'}/$gg{'pref'}_comb.bam
	my @fixBams; 
	for my $k1 (qw/tgt_fa src_fa/) {
		for (my $i=0; $i<@{$gg{$k1}}; $i++) {
			&fileSunhh::write2file(
				"$gg{'wrk_dir'}/$gg{'pref'}.pref_RS.${k1}_$i", 
				join("\t", 
					"${k1}_$i", 
					"${k1}_$i", 
					"${k1}_$i", 
					"$gg{'pref'}.${k1}_$i", 
					"$gg{'in_fq1'}", 
					"NA", 
					"illumina", 
					"NA", 
					"NA"
				)."\n", 
				'>'
			); 
			$cmd = "perl $gg{'pl_runHisat2_with2pass'} -in_pref_list $gg{'wrk_dir'}/$gg{'pref'}.pref_RS.${k1}_$i -db_hisat2 $gg{$k1}[$i] -wrk_dir $gg{'wrk_dir'} -para_hisat2 \" -p 30 --dta -q --phred33 \" -max_mismat_cnt $gg{'max_mismat_cnt'} -max_mismat_ratio $gg{'max_mismat_ratio'}"; 
			&exeCmd_1cmd($cmd) and &stopErr("[CMD] Failed at cmd: $cmd\n"); 

			chdir($gg{'abs_wrk_dir'}); 
			$cmd = "perl $gg{'pl_fix_NHnum'} -inBam $gg{'pref'}.${k1}_${i}_srt.bam -outBam $gg{'pref'}.${k1}_${i}_fixNH.bam -exe_samtools $gg{'exe_samtools'} -para_samTsrt \" -\@ 10 -m 4G \""; 
			&exeCmd_1cmd($cmd) and &stopErr("[CMD] Failed at cmd: $cmd\n"); 
			chdir($gg{'abs_cur_dir'}); 

			push(@fixBams, "$gg{'wrk_dir'}/$gg{'pref'}.${k1}_${i}_fixNH.bam"); 
			&fileSunhh::_rmtree("$gg{'wrk_dir'}/$gg{'pref'}.${k1}_${i}_srt.bam"); 
		}
	}
	$cmd = "$gg{'exe_samtools'} merge -n $gg{'outBam'} "; 
	$cmd .= join(" ", @fixBams); 
	&exeCmd_1cmd($cmd) and &stopErr("[CMD] Failed at cmd: $cmd\n"); 

	# If I get here, everything is good. So I want to remove temporary files. 
	for my $k1 (qw/tgt_fa src_fa/) {
		for (my $i=0; $i<@{$gg{$k1}}; $i++) {
			&fileSunhh::_rmtree("$gg{'wrk_dir'}/$gg{'pref'}.pref_RS.${k1}_$i"); 
			&fileSunhh::_rmtree("$gg{'wrk_dir'}/$gg{'pref'}.${k1}_${i}_fixNH.bam"); 
		}
	}

	return; 
}# step1_run_aln() 

sub setGlob {
	$gg{'pref'}         = 'out'; 
	$gg{'in_fq1'}       = 'NA'; 
	$gg{'in_fq2'}       = 'NA'; 

	$gg{'exe_samtools'} = 'samtools'; 
	$gg{'exe_hisat2'}   = 'hisat2'; 
	$gg{'pl_runHisat2_with2pass'} = '/home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/runHisat2_with2pass.pl'; 
	$gg{'pl_fix_NHnum'}           = '/home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/fix_NHnum.pl'; 
	$gg{'max_mismat_cnt'}   = -1; 
	$gg{'max_mismat_ratio'} = -1; 
	$gg{'wrk_dir'}          = './'; 

	$gg{'abs_cur_dir'}     = &fileSunhh::_abs_path("./"); 
	$gg{'abs_wrk_dir'}     = &fileSunhh::_abs_path($gg{'wrk_dir'}); 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# perl $0   -pref $gg{'pref'}   -src_fa fa   -tgt_fa fa   -in_fq1 in_RS.fq.gz
#
#   -pref       [$gg{'pref'}] Required. Output prefix 
#     -outBam     [$gg{'wrk_dir'}/$gg{'pref'}_comb.bam] Not necessary. This is the final output bam. 
#   -wrk_dir    [$gg{'wrk_dir'}] Directory storing result files. 
#   -src_fa     [faFile] \@ Indexed fasta files for source taxa where transmitted reads come from. 
#   -tgt_fa     [faFile] \@ Indexed fasta files for target taxa where transmitted reads go to. 
#   -in_fq1     [faFile] Sample fastq to be processed. 
#   
#   -max_mismat_cnt    [$gg{'max_mismat_cnt'}]   Parameter for pl_runHisat2_with2pass. -1 means no limitation. 
#   -max_mismat_ratio  [$gg{'max_mismat_ratio'}] Parameter for pl_runHisat2_with2pass. -1 means no limitation. 
#
#
#   -exe_samtools      [$gg{'exe_samtools'}] 
#   -exe_hisat2        [$gg{'exe_hisat2'}]
#   -pl_runHisat2_with2pass    [$gg{'pl_runHisat2_with2pass'}]
#   -pl_fix_NHnum              [$gg{'pl_fix_NHnum'}]
################################################################################
HH

	return; 
}# setGlob() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 
	for my $k (qw/in_fq1 src_fa tgt_fa/) {
		unless ( defined $opts{$k} ) {
			&tsmsg("[Err] -$k is required.\n\n"); 
			&LogInforSunhh::usage($gg{'help_txt'}); 
		}
	}

	# Outer tools
	for my $k (qw/exe_samtools exe_hisat2 pl_runHisat2_with2pass pl_fix_NHnum/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	for my $k (qw/pref in_fq1 max_mismat_cnt max_mismat_ratio wrk_dir src_fa tgt_fa/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	$gg{'outBam'}           = "$gg{'wrk_dir'}/$gg{'pref'}_comb.bam"; 
	defined $opts{'outBam'} and $gg{'outBam'} = $opts{'outBam'}; 

	$gg{'abs_wrk_dir'} = &fileSunhh::_abs_path($gg{'wrk_dir'}); 

	return; 
}# applyOpt() 

