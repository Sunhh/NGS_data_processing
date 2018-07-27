#!/usr/bin/perl -w
# 2018-07-27 : Extract sub-function as a perl script for easier use. 
# data processing steps for DNA sequencing I want to use : 

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
	"intervalLen:i",  # -1. If this is bigger than 0, I will combine GVCFs with interval list with multi-threads. 
); 

################################################################################
##########    Setting basic parameters 
################################################################################
# check if the input is sufficient to go on. 
my %gg; # global variates;  
my %cfg; 
&input_good() or &LogInforSunhh::usage($gg{'usage_txt'}); 
&set_pm(); # Setup multi-threading; 
&set_stepPara(); # Set parameter for detailed steps; 

### Method with CombineGVCFs : 
###   For most of gatk3 versions, I have to use CombineGVCFs instead of the following GenomicsDBImport to aggregate the GVCF files. 
###   Since the speed of combining is too slow, I want to use intervel lists to do it in parallel. 
$gg{'gvcf_list'} = [ map { $_->[0] } &fileSunhh::load_tabFile( $opts{'in_gvcf_list'} ) ]; 
$gg{'gvcf_num'}  = scalar( @{$gg{'gvcf_list'}} ); 
$cfg{'No_combineGVCF'} //= 50; 
&fileSunhh::write2file("$opts{'prj_ID'}_jnVCF.list", '', '>'); 

for (my $i=0; $i<$gg{'gvcf_num'}; $i+=$cfg{'No_combineGVCF'}) {
	$gg{'gvcf_jIdx'} ++; # Joined GVCF number; 
	my @curr_glist; 
	for (my $j=$i; $j<$gg{'gvcf_num'} and $j<$i+$cfg{'No_combineGVCF'}; $j++) {
		push(@curr_glist, $gg{'gvcf_list'}[$j]); 
	}
	my $jnGVCF ; 
	$jnGVCF = &step8_combineGVCF_interval( \@curr_glist, "$opts{'prj_ID'}_jn$gg{'gvcf_jIdx'}", '', "step8_cmd.combineGVCF.$opts{'prj_ID'}_$gg{'gvcf_jIdx'}", $gg{'wrk_dir'} ); 
	&fileSunhh::write2file("$opts{'prj_ID'}_jnVCF.list", "$jnGVCF\n", '>>'); 
}

### Step 8. Call variants in combined GVCFs. (rawV.vcf)
################################################################################
##########    Sub-routines 
################################################################################

sub input_good {
	$gg{'usage_txt'} = <<'H1'; 
################################################################################
perl $0 fasdfasf

  This program is used to run GATK-step8 : Generate combined_GVCF from input GVCF_list. 

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
  For CombineGVCFs
    -intervalLen        [number] Default -1. If this is bigger than 0, I will combine GVCFs with interval list with multi-threads. 

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

# Set multi-threads; 
sub set_pm {
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

# Set parameter for detailed steps; 
sub set_stepPara {
	# Set parameter for detailed steps; 
	# Get the path for diction of ref_fasta; 
	$cfg{'ref_fasta'} = &fileSunhh::_abs_path_4link( $cfg{'ref_fasta'} ); # Use the full path of 'ref_fasta'; 
	&generate_ref_idx( $cfg{'ref_fasta'} ) or &stopErr("[Err] Cannot generate enough index files for ref_fasta: $cfg{'ref_fasta'}\n"); 
	# Set interval length : 
	$gg{'para'}{'intervalLen'} = -1; 
	defined $cfg{'intervalLen'} and $gg{'para'}{'intervalLen'} = $cfg{'intervalLen'}; 
	defined $opts{'intervalLen'} and $gg{'para'}{'intervalLen'} = $opts{'intervalLen'}; 
}# set_stepPara() 


sub step8_combineGVCF_interval {
	# Return the combined GVCF file's name. 
	# Using interval list to combine GVCF in parallel, and then use CombineVariants to join them. 
	#   The speed of 'CombineVariants' is also slow, so I want to use my own perl script to do the combination, and compress it with tabix. 
	my ($glist, $opref, $interval_list, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my $cmd = ''; 
	if ( $interval_list ne '' or $gg{'para'}{'intervalLen'} <= 0 ) {
		$cmd .= "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_gatk'} -T CombineGVCFs \\\n"; # Add command here; 
		$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
		$interval_list ne '' and $cmd .= "  -L $interval_list \\\n"; 
		$cmd .= "  -o ${opref}.g.vcf.gz \\\n"; 
		for my $gFn (@$glist) {
			$cmd .= "  -V $gFn \\\n"; 
		}
		$cmd .= "1> s8.std.combineGVCF.${opref} 2> s8.err.combineGVCF.${opref}\n"; 
	} else {
		# Get interval lists; 
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
		my @subPref_list; 
		for (my $i=0; $i<@seqLen; $i++) {
			my ($id, $len) = @{$seqLen[$i]};
			for (my $j=1; ($j-1)*$gg{'para'}{'intervalLen'}+1 < $len; $j++) {
				my $sub_gvcf_pref = "${opref}.interval_result.$i.$j"; 
				push(@subPref_list, $sub_gvcf_pref); 
				my $s = ($j-1)*$gg{'para'}{'intervalLen'}+1; 
				my $e = $s + $gg{'para'}{'intervalLen'}-1; 
				$e > $len and $e = $len; 
				my $ti_fn = "${opref}.$i.$j.interval_list"; 

				$gg{'MAX_PROCESSES'} = &LogInforSunhh::change_procN( $gg{'pm'}, "$gg{'nprocF'}", $gg{'MAX_PROCESSES'} ); 
				my $pid = $gg{'pm'}->start and next; 
				&fileSunhh::write2file($ti_fn, "$header_txt" . join("\t", $id, $s, $e, '+', "${id}_${s}_${e}"). "\n", '>'); # It doesn't matter if strand is '+'/'-' for CombineGVCF; 
				&step8_combineGVCF_interval( $glist, $sub_gvcf_pref, $ti_fn, "${fn_cmd}.interval_cmd.$i.$j", "./" ); 
				$gg{'pm'}->finish; 
			}
		}
		$gg{'pm'}->wait_all_children; 
		if ( $opts{'plCatVar'} ) {
			my $subGVCF_lisFn = "list_subGVCF.${opref}"; 
			&fileSunhh::write2file($subGVCF_lisFn, '', '>'); 
			for my $tfP (@subPref_list) {
				&fileSunhh::write2file($subGVCF_lisFn, "${tfP}.g.vcf.gz\n", '>>'); 
			}
			$cmd .= "$cfg{'exe_perl'} $cfg{'pl_catVar'} $subGVCF_lisFn | \\\n"; 
			$cmd .= "  $cfg{'exe_bgzip'} -c /dev/stdin \\\n"; 
			$cmd .= "1> ${opref}.g.vcf.gz\n"; 

			$cmd .= "$cfg{'exe_tabix'} -p vcf ${opref}.g.vcf.gz\n"; 
		} else {
			$cmd .= "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T CombineVariants \\\n"; 
			$cmd .= "  -nt 20 \\\n"; 
			$cmd .= "  --assumeIdenticalSamples -genotypeMergeOptions UNSORTED \\\n"; 
			$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
			$cmd .= "  -o ${opref}.g.vcf.gz \\\n"; 
			for my $tfP (@subPref_list) {
				$cmd .= "  -V ${tfP}.g.vcf.gz \\\n"; 
			}
			$cmd .= "1> s8.std.combineGVCF.${opref} 2> s8.err.combineGVCF.${opref}\n"; 
		}
		for my $tfP (@subPref_list) {
			$cmd .= "rm -f ${tfP}.g.vcf.gz ${tfP}.g.vcf.gz.tbi\n"; 
		}
	}
	
	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 8 commands: Combine GVCFs.\n" ); 
		chdir($ori_dir); 
		return( "${opref}.g.vcf.gz" ); 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return("${opref}.g.vcf.gz"); 
}# step8_combineGVCF_interval() 


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

# Generate reference index; 
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

