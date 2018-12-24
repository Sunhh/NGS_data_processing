#!/usr/bin/perl -w
# 2018-08-01 : Drop use of ' -nt '
# 2018-07-28 : Run GATK step9 separately. 
### Step 9. call_rawV + filter_rawSNP + filter_rawInDel + Combine_filtered_SNP_and_InDel + slct_passed_data; This can be done for each CHR. 
#
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use ConfigSunhh; 
my $cfgs_obj = ConfigSunhh->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"in_gvcf_list:s", # Required. Format : GVCF_file_1\nGVCF_file_2\n
	"prj_ID:s",       # 'outGATK'. A project ID for merging GVCF files. 
	"out_dir:s",      # Output dir; 
	"wrk_dir:s",      # Work dir; 
	
	"conf_file:s",    # Required. This file tells the path information of softwares. 
	"cpuN:i",         # 
	"plCatVar!",      # 

	# Detail parameters 
	"CallByScf!",     # 
	"intervalLen:i",  # -1. Used in step8, and I want to use this in step9. 
); 

################################################################################
##########    Setting basic parameters 
################################################################################
# check if the input is sufficient to go on. 
my %gg; # global variates;  
my %cfg; 
&input_good() or &LogInforSunhh::usage($gg{'usage_txt'}); 
&set_pm(); # Set multi-threads; 
&set_stepPara(); # Set parameter for detailed steps; 

# Section three : Call variants and filter data for the first time; 
### It should be possible to do the steps one by one for each scfID : 
###   call_rawV + filter_rawSNP + filter_rawInDel + Combine_filtered_SNP_and_InDel + slct_passed_data; 
&step9_gvcf2var( $gg{'jnGVCF_list'}, "$opts{'prj_ID'}", '', "step9_cmd.callRawV.$opts{'prj_ID'}", $gg{'wrk_dir'} ); 


################################################################################
##########    Sub-routines 
################################################################################

sub input_good {
	$gg{'usage_txt'} = <<'H1'; 
################################################################################
perl $0 fasdfasf

  This program is used to run step9_gvcf2var in pipe_gatk separately. 

  -help                 Show this help. 

  -conf_file            [filename] Required. Such as 'pipe_gatk_conf'; 
                          This file tells the path information of softwares.
  -in_gvcf_list         [filename] Required. Format : GVCF_file_1\\nGVCF_file_2\\n
  -prj_ID               [String] Output prefix (default 'outGATK') for the whole project to merge GVCF files. Cannot contain directory; 

  -wrk_dir              [dirname] Assign temporary directory to work in. 
  -cpuN                 [1]       Number of samples to be processed in parallel. 
  -plCatVar             [Boolean] If given, I will use perl instead of 'CombineVariants' to join the interval-called GVCFs. 
                                  The 'CatVariants' is not available in my GATK version. 

  ### Specific parameters : These parameters will over-write conf_file. 
  For GenotypeGVCFs     
    -CallByScf          [Boolean] If given, I'll try to call variants for each scaffold with multi-threads; 
    -intervalLen        [number] Default -1. If this is bigger than 0, I will run genotyping with intervals when not breaking scaffolds. 

################################################################################
H1

	defined $opts{'conf_file'} or return 0; 
	defined $opts{'in_gvcf_list'} or return 0; 
	$opts{'prj_ID'} //= 'outGATK'; 
################################################################################
##########    Example configure file
################################################################################
$gg{'example'}{'text_cfg'} = <<'CFG'; 

pl_getSam             /home/Sunhh/tools/github/NGS_data_processing/Proc_Sam/get_required_sam.pl
pl_batchRun           /home/Sunhh/tools/github/NGS_data_processing/run_cmd_in_batch.pl
batchNum              1

pl_catVar             /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/gatk/CatVariants.pl
exe_bgzip             /opt/align/samtools/install/htslib-1.5/bin/bgzip
exe_tabix             /opt/align/samtools/install/htslib-1.5/bin/tabix

exe_perl              /usr/bin/perl
exe_java              /usr/java/jre1.8.0_144/bin/java
#exe_java              /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.141-1.b16.fc26.x86_64/jre/bin/java
#JAVA_HOME             /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.141-1.b16.fc26.x86_64/
jar_picard            /home/Sunhh/src/align/picard/v2.10.3/picard.jar
jar_gatk              /Data/Sunhh/src/pipeline/gatk/nightly_20171018/GenomeAnalysisTK.jar
exe_samtools          /opt/align/samtools/samtools-1.5/samtools
para_maoPU            -q 17 -Q 0 --ff 3844
# flag 3844 includes 'unmapped', 'not primary alignment/secondary alignment', 'read fails platform/vendor quality checks', 'read is PCR or optical duplicate', 'supplementary alignment'
exe_bwa               /opt/align/bwa/bwa-0.7.15/bwa
para_bwa              -M -t 70

ref_fasta             /Data/Sunhh/wm_reseqByPB/db/wm97pb_v2ID.scf.fa

dir_tmp               /Data/Sunhh/wm_reseqByPB/01_get_gvcf_pbv2/temp

nct_hapCaller         5


dirMao                /home/Sunhh/tools/github/NGS_data_processing/reseq_tools/mao_exe
exe_combine2PileFiles __dirMao__/combine2PileFiles
exe_reSeqPrintRefChr  __dirMao__/reSeqPrintRefChr
exe_reSeqPrintSample  __dirMao__/reSeqPrintSample.indel.fast.strAssign
para_reSeqPrintSample 2 2 0.3
#exe_reSeqPrintSample  __dirMao__/reSeqPrintSample.indel.fast
#exe_reSeqPrintSample  __dirMao__/reSeqPrintSample.indel.fast.strAssign.moreHeter
#exe_rmRedunSam2       __dirMao__/rmRedunSam2
#exe_rmRedunSam3       __dirMao__/rmRedunSam3
No_combineGVCF         50

exe_bcftools           /opt/align/samtools/bcftools-1.5/bcftools
grp_sbSampleN          40

dirDbsnp              /data/Sunhh/watermelon/01.source_reseq/01.gatk/independent_snps
known_dbsnp            -knownSites __dirDbsnp__/WmReseq_sample464_rmClose.sMao_dbsnp.vcf.gz  -knownSites __dirDbsnp__/WmReseq_sample464.sSB_dbsnp.vcf.gz 
para_bqsr              -nct 40

CFG

################################################################################
##########    Example configure file
################################################################################

	# Setup directories to be used. 
	if ( defined $opts{'wrk_dir'} ) {
		$gg{'wrk_dir'} = &fileSunhh::_abs_path($opts{'wrk_dir'}); 
		-e $gg{'wrk_dir'} or do { mkdir($gg{'wrk_dir'}) or &stopErr("[Err] Failed to create wrk_dir [$gg{'wrk_dir'}]\n"); }; 
	} else {
		$gg{'wrk_dir'} = &fileSunhh::new_tmp_dir('create' => 1) or &stopErr("[Err] Failed to create tmp dir\n"); 
		$gg{'wrk_dir'} = &fileSunhh::_abs_path($gg{'wrk_dir'}); 
	}
	$gg{'ori_dir'} = &fileSunhh::_abs_path("./"); 
	$opts{'out_dir'} //= 'gatk_out'; 
	$gg{'out_dir'} = $opts{'out_dir'}; 

	# Setup configuration parameters. 
	&fileSunhh::write2file("$gg{'wrk_dir'}/example.conf", $gg{'example'}{'text_cfg'}, '>'); 
	$cfgs_obj->getConfig('cfg_file'=>"$gg{'wrk_dir'}/example.conf", 'replace'=>0, 'hash_r'=>\%cfg); 
	$cfgs_obj->getConfig('cfg_file'=>$opts{'conf_file'}, 'replace'=>1, 'hash_r'=>\%cfg); 
	$cfgs_obj->writeConfig('cfg_file'=>"$gg{'wrk_dir'}/$opts{'prj_ID'}_gatk.conf", 'hash_r'=>\%cfg); 
	$cfg{'dir_tmp'} = &fileSunhh::_abs_path($cfg{'dir_tmp'}); 
	unless ( -e $cfg{'dir_tmp'} ) {
		mkdir($cfg{'dir_tmp'}) or &stopErr("[Err] Failed to create dir_tmp [$cfg{'dir_tmp'}]\n"); 
	}
	return 1; 
}# input_good() 

sub set_pm {
	# Set multi-threads; 
	$gg{'MAX_PROCESSES'} = ( defined $opts{'cpuN'} ) 
		? $opts{'cpuN'} 
		: ( defined $cfg{'batchNum'} ) 
			? $cfg{'batchNum'} 
			: 1
	; 
	$gg{'nprocF'} = "$opts{'prj_ID'}.Nproc"; # This file will exist in wrk_dir; 
	&fileSunhh::write2file("$gg{'wrk_dir'}/$gg{'nprocF'}", "$gg{'MAX_PROCESSES'}\n", '>'); 
	$gg{'pm'} = &LogInforSunhh::get_pm( $gg{'MAX_PROCESSES'} ); 
}# set_pm() 

sub set_stepPara {
	# Set parameter for detailed steps; 
	# Get the path for diction of ref_fasta; 
	$cfg{'ref_fasta'} = &fileSunhh::_abs_path_4link( $cfg{'ref_fasta'} ); # Use the full path of 'ref_fasta'; 
	&generate_ref_idx( $cfg{'ref_fasta'} ) or &stopErr("[Err] Cannot generate enough index files for ref_fasta: $cfg{'ref_fasta'}\n"); 

	# Set interval length : 
	$gg{'para'}{'CallByScf'} = 0; 
	defined $cfg{'CallByScf'} and $gg{'para'}{'CallByScf'} = $cfg{'CallByScf'}; 
	$opts{'CallByScf'} and $gg{'para'}{'CallByScf'} = 1; 
	$gg{'para'}{'CallByScf'} =~ s!^\s*(false|F)\s*$!0!i; # True/False/0/1; 
	$gg{'para'}{'CallByScf'} =~ s!^\s*(true|T)\s*$!1!i; # True/False/0/1; 
	$gg{'para'}{'intervalLen'} = -1; 
	defined $cfg{'intervalLen'}  and $gg{'para'}{'intervalLen'} = $cfg{'intervalLen'}; 
	defined $opts{'intervalLen'} and $gg{'para'}{'intervalLen'} = $opts{'intervalLen'}; 

	# Get GVCF list from in_gvcf_list; 
	$gg{'jnGVCF_list'} = [ map { $_->[0] } &fileSunhh::load_tabFile( $opts{'in_gvcf_list'} ) ]; 

}# set_stepPara() 

sub generate_ref_idx {
	my $fn = shift;
	my $pref = $fn;
	$pref =~ s!\.fasta$|\.fa$|\.fas$!! or &stopErr("[Err] ref_fasta [$fn] should be ended with .fasta/.fa\n");
	unless (-e "${fn}.fai") {
		&exeCmd_1cmd("$cfg{'exe_samtools'} faidx $fn") and return 0; # here 'return 0' means this program failed to run correctly.
		-e "${fn}.fai" or return 0;
	}
	unless (-e "${pref}.dict") {
		&exeCmd_1cmd("$cfg{'exe_java'} -jar $cfg{'jar_picard'} CreateSequenceDictionary   R=$fn   O=${pref}.dict") and return 0;
		-e "${pref}.dict" or return 0;
	}
	$gg{'ref_fai'}  = "${fn}.fai"; 
	$gg{'ref_dict'} = "${pref}.dict"; 
	return 1;
}# generate_ref_idx()

sub step9_gvcf2var {
	# call_rawV + filter_rawSNP + filter_rawInDel + Combine_filtered_SNP_and_InDel + slct_passed_data; 
	# I don't know if it is OK to set the interval boundary overlapping an InDel variant, so I want to do this calling within any whole scaffold. 
	# Ref : https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
	# Ref : https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
	# I need to check this ref : https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/JointGenotypingWf.wdl
	
	my ($glist, $opref, $scfID, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my $cmd = ''; 

	if ( $scfID ne '' or !$gg{'para'}{'CallByScf'} ) {
		$cmd .= "$cfg{'exe_java'} -Xmx10g -Xms10g -jar $cfg{'jar_gatk'} -T GenotypeGVCFs \\\n"; 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$cmd .= "  -o ${opref}_rawV.vcf.gz \\\n"; 
		$scfID ne '' and $cmd .= "  -L $scfID \\\n"; 
		for my $gFn (@$glist) {
			$cmd .= "  --variant $gFn \\\n"; 
		}
		$cmd .= "1> s9.std.callRawV.${opref} 2> s9.err.callRawV.${opref}\n"; 

		$cmd .= "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T SelectVariants \\\n"; 
		$cmd .= "  -selectType SNP \\\n"; 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$cmd .= "  -V ${opref}_rawV.vcf.gz \\\n"; 
		$cmd .= "  -o ${opref}_rawV_snps.vcf.gz \\\n"; 
		$cmd .= "1> s9.std.getRawSNP.${opref} 2> s9.err.getRawSNP.${opref}\n"; 

		$cmd .= "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T VariantFiltration \\\n"; 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$cmd .= "  -V ${opref}_rawV_snps.vcf.gz \\\n"; 
		$cmd .= "  -o ${opref}_filtV_snps.vcf.gz \\\n"; 
		$cmd .= "  --filterExpression \"QD < 2.0\"              --filterName \"QD_snp\" \\\n";
		$cmd .= "  --filterExpression \"FS > 60.0\"             --filterName \"FS_snp\" \\\n";
		$cmd .= "  --filterExpression \"MQ < 40.0\"             --filterName \"MQ_snp\" \\\n";
		$cmd .= "  --filterExpression \"MQRankSum < -12.5\"     --filterName \"MQRankSum_snp\" \\\n";
		$cmd .= "  --filterExpression \"ReadPosRankSum < -8.0\" --filterName \"ReadPosRankSum_snp\" \\\n";
		$cmd .= "1> s9.std.getFiltSNP.${opref} 2> s9.err.getFiltSNP.${opref}\n"; 

		$cmd .= "$cfg{'exe_java'} -Xmx10G  -jar $cfg{'jar_gatk'} -T SelectVariants \\\n"; 
		$cmd .= "  -selectType INDEL \\\n"; 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$cmd .= "  -V ${opref}_rawV.vcf.gz \\\n"; 
		$cmd .= "  -o ${opref}_rawV_indels.vcf.gz \\\n"; 
		$cmd .= "1> s9.std.getRawInDel.${opref} 2> s9.err.getRawInDel.${opref}\n"; 

		$cmd .= "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T VariantFiltration \\\n"; 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$cmd .= "  -V ${opref}_rawV_indels.vcf.gz \\\n"; 
		$cmd .= "  -o ${opref}_filtV_indels.vcf.gz \\\n"; 
		$cmd .= "  --filterExpression \"QD < 2.0\"               --filterName \"QD_indel\" \\\n";
		$cmd .= "  --filterExpression \"FS > 200.0\"             --filterName \"FS_indel\" \\\n";
		$cmd .= "  --filterExpression \"ReadPosRankSum < -20.0\" --filterName \"ReadPosRankSum_indel\" \\\n";
		$cmd .= "1>s9.std.getFiltInDel.${opref} 2> s9.err.getFiltInDel.${opref}\n"; 

		$cmd .= "$cfg{'exe_java'} \\\n"; 
		$cmd .= "  -Djava.io.tmpdir=$cfg{'dir_tmp'} \\\n"; 
		$cmd .= "  -jar $cfg{'jar_gatk'} -T CombineVariants \\\n"; 
		# $cmd .= "  -nt 4 \\\n"; # I drop this parameter because it causes failure sometimes. 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$cmd .= "  -o ${opref}_filtV.vcf.gz \\\n"; 
		$cmd .= "  --variant:snp    ${opref}_filtV_snps.vcf.gz \\\n"; 
		$cmd .= "  --variant:indel  ${opref}_filtV_indels.vcf.gz \\\n"; 
		$cmd .= "  -genotypeMergeOptions PRIORITIZE   -priority snp,indel \\\n"; 
		$cmd .= "1>s9.std.combineFiltV.${opref} 2> s9.err.combineFiltV.${opref}\n"; 

	} else {
		my $ifh = &openFH($gg{'ref_dict'}, '<'); 
		my (@seqLen, $header_txt); 
		while (<$ifh>) {
			chomp; 
			$header_txt .= "$_\n"; 
			m!^\@SQ\t! or next; 
			my @ta = split(/\t/, $_); 
			$ta[1] =~ m!^SN:(\S+)$! or &stopErr("[Err] Unknown ID [$ta[1]]\n"); 
			push(@seqLen, [$1]); 
			$ta[2] =~ m!^LN:(\d+)$! or &stopErr("[Err] Unknown len [$ta[2]]\n"); 
			$seqLen[-1][1] = $1; 
		}
		close($ifh); 
		@seqLen = sort { $b->[1] <=> $a->[1] } @seqLen; 
		my @subPref_list; 
		my @toRun_intervals; # ( [interval_scfIDs, interval_total_length], [], ... )
		for (my $i=0; $i<@seqLen; $i++) {
			my ($id, $len) = @{$seqLen[$i]}; 
			if (@toRun_intervals == 0) {
				push(@toRun_intervals, [ [[$id, $len]], $len ]); 
			} elsif ( $toRun_intervals[-1][1] >= $gg{'para'}{'intervalLen'} ) {
				push(@toRun_intervals, [ [[$id, $len]], $len ]); 
			} else {
				push(@{$toRun_intervals[-1][0]}, [$id, $len]); 
				$toRun_intervals[-1][1] += $len; 
			}
		}
		&tsmsg("[Rec] Total " . scalar(@toRun_intervals) . " intervals will be processed parallelly.\n"); 
		for (my $i=0; $i<@toRun_intervals; $i++) {
			my $sub_gvcf_pref = "${opref}.gvcf2var_sepChr.$i"; 
			push(@subPref_list, $sub_gvcf_pref); 
			my $ti_fn = "${sub_gvcf_pref}.interval_list"; 
			&fileSunhh::write2file($ti_fn, "$header_txt", '>'); 
			for my $iL (@{$toRun_intervals[$i][0]}) {
				&fileSunhh::write2file($ti_fn, join("\t", $iL->[0], 1, $iL->[1], "+", join("_", $iL->[0], 1, $iL->[1]))."\n", '>>'); 
			}
			$gg{'MAX_PROCESSES'} = &LogInforSunhh::change_procN( $gg{'pm'}, "$gg{'nprocF'}", $gg{'MAX_PROCESSES'} ); 
			my $pid = $gg{'pm'}->start and next; 
			&step9_gvcf2var( $glist, $sub_gvcf_pref, $ti_fn, "${fn_cmd}.gvcf2var_cmd.$i", "./" ); 
			$gg{'pm'}->finish; 
		}
		$gg{'pm'}->wait_all_children; 
		if ( $opts{'plCatVar'} ) {
			my $subGVCF_lisFn1 = "list_subFiltV.${opref}"; 
			my $subGVCF_lisFn2 = "list_subrawV.${opref}"; 
			&fileSunhh::write2file($subGVCF_lisFn1, '', '>'); 
			&fileSunhh::write2file($subGVCF_lisFn2, '', '>'); 
			for my $tfP (@subPref_list) {
				&fileSunhh::write2file($subGVCF_lisFn1, "${tfP}_filtV.vcf.gz\n", '>>'); 
				&fileSunhh::write2file($subGVCF_lisFn2, "${tfP}_rawV.vcf.gz\n", '>>'); 
			}
			$cmd .= "$cfg{'exe_perl'} $cfg{'pl_catVar'} $subGVCF_lisFn1 | \\\n"; 
			$cmd .= "  $cfg{'exe_bgzip'} -c /dev/stdin \\\n"; 
			$cmd .= "1> ${opref}_filtV.vcf.gz\n"; 

			$cmd .= "$cfg{'exe_tabix'} -p vcf ${opref}_filtV.vcf.gz\n"; 

			$cmd .= "$cfg{'exe_perl'} $cfg{'pl_catVar'} $subGVCF_lisFn2 | \\\n"; 
			$cmd .= "  $cfg{'exe_bgzip'} -c /dev/stdin \\\n"; 
			$cmd .= "1> ${opref}_rawV.vcf.gz\n"; 

			$cmd .= "$cfg{'exe_tabix'} -p vcf ${opref}_rawV.vcf.gz\n"; 
		} else {
			$cmd .= "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T CombineVariants \\\n"; 
			# $cmd .= "  -nt 20 \\\n"; 
			$cmd .= "  --assumeIdenticalSamples -genotypeMergeOptions UNSORTED \\\n"; 
			$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
			$cmd .= "  -o ${opref}_filtV.vcf.gz \\\n"; 
			for my $tfP (@subPref_list) {
				$cmd .= "  -V ${tfP}_filtV.vcf.gz \\\n"; 
			}
			$cmd .= "1> s9.std.gvcf2var_filtV.${opref} 2> s9.err.gvcf2var_filtV.${opref}\n"; 

			$cmd .= "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T CombineVariants \\\n"; 
			# $cmd .= "  -nt 20 \\\n"; 
			$cmd .= "  --assumeIdenticalSamples -genotypeMergeOptions UNSORTED \\\n"; 
			$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
			$cmd .= "  -o ${opref}_rawV.vcf.gz \\\n"; 
			for my $tfP (@subPref_list) {
				$cmd .= "  -V ${tfP}_rawV.vcf.gz \\\n"; 
			}
			$cmd .= "1> s9.std.gvcf2var_rawV.${opref} 2> s9.err.gvcf2var_rawV.${opref}\n"; 
		}
		for my $tfP (@subPref_list) {
			$cmd .= "rm -f ${tfP}.vcf.gz ${tfP}.vcf.gz.tbi\n"; 
		}
	}

	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 2 commands: Mark illumina adapters\n" ); 
		chdir($ori_dir); 
		return("${opref}.vcf.gz"); 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return("${opref}.vcf.gz"); 
}# step9_gvcf2var() 


sub stdout_file {
	my ($fn, @txt) = @_; 
	for (@txt) {
		print STDOUT $_; 
	}
	my $ifh = &openFH($fn, '<'); 
	while (<$ifh>) {
		print STDOUT $_; 
	}
	close($ifh); 
}# stdout_file() 

# Group @fq_infor by samples; 
# {RGSM}{'fq_infor'} => [ $fq_infor[$i], $fq_infor[$i], ... ]; 
# {RGSM}{'in_order'} => number; 
sub groupFqBySM {
	my ($fqR) = @_; 
	my @fq_infor = @$fqR; 
	my %sm_infor; 
	for (my $i=0; $i<@fq_infor; $i++) {
		my %fqHash = %{$fq_infor[$i]}; 
		push(@{$sm_infor{$fqHash{'SM'}}{'fq_infor'}}, $fq_infor[$i]); 
		$sm_infor{$fqHash{'SM'}}{'in_order'} //= $i; 
	}
	return(\%sm_infor); 
}# groupFqBySM() 

# sub tt {
# 	my ($smID, $fn_cmd, $wrk_dir) = @_; 
# 	my $ori_dir = &fileSunhh::_abs_path("./"); 
# 	chdir($wrk_dir); 
# 	my $cmd = ''; 
# 	$cmd .= ""; # Add command here; 
#
# 	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
# 	if ($gg{'cmd'}) {
# 		&stdout_file( $fn_cmd, "# Step 2 commands: Mark illumina adapters\n" ); 
# 		chdir($ori_dir); 
# 		return; 
# 	}
# 	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
# 	chdir($ori_dir); 
# 	return; 
# }# tt() 

