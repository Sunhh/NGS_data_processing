#!/usr/bin/perl
# 2018-10-25 : step X1 . This is used to align non/self-grafted rnaseq samples to grafted-source genome, and get read alignments for later filtering of possible transmitted reads. 
#   If a possible transmitted read is included by a aligned non/self-grafted read, this is a false positive one. 
#   Input  : -src_fa , -in_fq1 :@
#   Output : pref_comb.nonTrans.bam
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
	"in_fq1:s@",  # fastq_R1; These are non/self-grafted samples. 
	"wrk_dir:s", # './', Assign a directory to store all resulting bam files, instead of in current folder.
	"outBam:s",  # The combined bam file output. This will replace the output file "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam" 

	"exe_samtools:s", 
	"exe_hisat2:s", 
	"pl_runHisat2_with2pass:s", 
		"max_mismat_cnt:i", 
		"max_mismat_ratio:f", 
); 

my %gg; 
&setGlob(); 
&applyOpt(); 
&stepX1_aln_tgtRd(); 

sub stepX1_aln_tgtRd {
	# Step X1 : Align non/self-grafted (target) samples to source genome, to get background for false-positive check. 
	# Produce : $gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam
	my $cmd = ''; 
	my @toCombBam; 
	for (my $j=0; $j<@{$gg{'in_fq1'}}; $j++) {
		for (my $i=0; $i<@{$gg{'src_fa'}}; $i++) {
			&fileSunhh::write2file(
				"$gg{'wrk_dir'}/$gg{'pref'}.pref_RS.$j.src_$i", 
				join("\t", 
					"src_${j}_$i", 
					"src_${j}_$i", 
					"src_${j}_$i", 
					"$gg{'pref'}.src_${j}_$i", 
					"$gg{'in_fq1'}[$j]", 
					"NA", 
					"illumina", 
					"NA", 
					"NA"
				), 
				'>'
			); 
			$cmd = "perl $gg{'pl_runHisat2_with2pass'} -in_pref_list $gg{'wrk_dir'}/$gg{'pref'}.pref_RS.$j.src_$i -db_hisat2 $gg{'src_fa'}[$i] -wrk_dir $gg{'wrk_dir'} -para_hisat2 \" -p 30 --dta -q --phred33 \" -max_mismat_cnt $gg{'max_mismat_cnt'} -max_mismat_ratio $gg{'max_mismat_ratio'}"; 
			&exeCmd_1cmd($cmd) and &stopErr("[CMD] Failed at cmd: $cmd\n"); 
			push(@toCombBam, "$gg{'wrk_dir'}/$gg{'pref'}.src_${j}_${i}_srt.bam"); 
		}
	}
	if (scalar(@toCombBam) > 1) {
		$cmd = "$gg{'exe_samtools'} merge $gg{'outBam'} "; 
		$cmd .= join(" ", @toCombBam); 
		&exeCmd_1cmd($cmd) and &stopErr("[CMD] Failed at cmd: $cmd\n"); 
	} elsif (scalar(@toCombBam) == 1) {
		&fileSunhh::_move( $toCombBam[0], $gg{'outBam'} ); 
	} else {
		&stopErr("[Err] No bam files generated?\n"); 
	}

	# If I get here, everything is good. So I want to remove temporary files. 
	for (my $j=0; $j<@{$gg{'in_fq1'}}; $j++) {
		for (my $i=0; $i<@{$gg{'src_fa'}}; $i++) {
			&fileSunhh::_rmtree("$gg{'wrk_dir'}/$gg{'pref'}.pref_RS.$j.src_$i"); 
			&fileSunhh::_rmtree("$gg{'wrk_dir'}/$gg{'pref'}.src_${j}_${i}_srt.bam"); 
		}
	}

	return; 
}# stepX1_aln_tgtRd() 

sub applyOpt {
	$opts{'help'} and &LogInforSunhh::usage($gg{'help_txt'}); 

	for my $k (qw/in_fq1 src_fa/) {
		unless ( defined $opts{$k} ) {
			&tsmsg("[Err] -$k is required.\n\n"); 
			&LogInforSunhh::usage($gg{'help_txt'}); 
		}
	}

	# For outer tools
	for my $k (qw/exe_samtools exe_hisat2 pl_runHisat2_with2pass/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	for my $k (qw/pref in_fq1 wrk_dir src_fa/) {
		defined $opts{$k} and $gg{$k} = $opts{$k}; 
	}

	$gg{'outBam'} = "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam"; 
	defined $opts{'outBam'} and $gg{'outBam'} = $opts{'outBam'}; 

	return; 
}# applyOpt() 

sub setGlob {
	$gg{'pref'}    = 'pref'; 
	$gg{'in_fq1'}  = []; 
	$gg{'wrk_dir'} = "./"; 
	$gg{'outBam'}  = "$gg{'wrk_dir'}/$gg{'pref'}_comb.nonTrans.bam"; 

	$gg{'max_mismat_cnt'}    = -1; 
	$gg{'max_mismat_ratio'}  = -1; 

	$gg{'exe_samtools'} = 'samtools'; 
	$gg{'exe_hisat2'}   = 'hisat2'; 
	$gg{'pl_runHisat2_with2pass'} = '/home/Sunhh/tools/github/NGS_data_processing/rnaseq_tools/runHisat2_with2pass.pl'; 

$gg{'help_txt'} = <<"HH"; 
################################################################################
# perl $0   -pref $gg{'pref'}   -src_fa fa   -in_fq1 in_RS1.fq.gz   [-in_fq1 in_RS2.fq.gz]
#
#   -pref       [$gg{'pref'}] Output prefix 
#     -outBam   [$gg{'outBam'}] Not necessary. 
#   -wrk_dir    [$gg{'wrk_dir'}] Directory storing result files. 
#   -src_fa     [faFile] \@ Indexed fasta files for source taxa where transmitted reads come from. 
#   -in_fq1     [faFile] \@ non/self-grafted sample fastq to be processed.
#
#   -max_mismat_cnt    [$gg{'max_mismat_cnt'}]   Parameter for pl_runHisat2_with2pass. -1 means no limitation.
#   -max_mismat_ratio  [$gg{'max_mismat_ratio'}] Parameter for pl_runHisat2_with2pass. -1 means no limitation.
#
#   -exe_samtools      [$gg{'exe_samtools'}]
#   -exe_hisat2        [$gg{'exe_hisat2'}]
#   -pl_runHisat2_with2pass    [$gg{'pl_runHisat2_with2pass'}]
################################################################################
HH
	return; 
}# setGlob() 


