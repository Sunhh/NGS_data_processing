#!/usr/bin/perl -w
# 2018-08-01 : I don't want to use ' -nt ' parameter any more! This often causes strange I/O errors! 
# 2018-07-05 : Update for help; 
# 2016-02-12 : Pipe GATK step by step for multiple input fastq files; 
# data processing steps for DNA sequencing I want to use : 
### Step 1. Change input fastq to raw 
### Step 2. Mark illumina adapters with picard MarkIlluminaAdapters; I change 'XT:' to 'YT:' for bwa information, because 'XT' is used in bwa; 
### Step 3. Align _mrkAdp.bam to reference fasta with BWA-MEM and merge it with uBam as aligned_pipe1.bam: out file "$gg{'wrk_dir'}/$fqHash{'pref'}_aln_pipe1.bam"; These files have all the reads and read base, and YT is recorded as adapters; 
### Step 4. Merge RGs into single sample file. 
### Step 5. Mark duplicates in bam file. out file : ${smID}_bySM_dedup.bam and index; 
### Step 6. This is reserved for further processing including BQSR , realignment, fix SetNmMdAndUqTags; out file : ${smID}_bySM_fix.bam ;
### Step 7. Call GVCF files for each sample using 'HaplotypeCaller'. It takes about 2.5-3 hours to process a bam to GVCF, so I don't bother to group CHRs for speed. 
### Step 8. Combine GVCFs. 
### Step 9. call_rawV + filter_rawSNP + filter_rawInDel + Combine_filtered_SNP_and_InDel + slct_passed_data; This can be done for each CHR. 
#

##### The following have not been applied. 
### Step 10. Do VQSR with both filtered SNP and indel on rawV.vcf . 
### Step 12. .... 
# 2018-06-19 Try to run all the processes in one command; 

use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use ConfigSunhh; 
my $cfgs_obj = ConfigSunhh->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"in_pref_list:s", # Required. Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2
	"prj_ID:s",       # 'outGATK'. A project ID for merging GVCF files. 
	"out_dir:s",      # Output dir; 
	"wrk_dir:s",      # Work dir; 
	
	"conf_file:s",    # Required. This file tells the path information of softwares. 
	"doStep:s", 
	"cpuN:i",         # 
	"plCatVar!",      # 

	# Detail parameters 
	"ERC:s",          # 
	"intervalLen:i",  # -1. If this is bigger than 0, I will combine GVCFs with interval list with multi-threads. 
	"CallByScf!",     # 
	"help!", 
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

# Section one   : Process data from fastq to GVCF by grouped sample; per-sample; 
for my $sm (sort { $gg{'smFq'}{$a}{'in_order'} <=> $gg{'smFq'}{$b}{'in_order'} } keys %{$gg{'smFq'}}) {
	$gg{'MAX_PROCESSES'} = &LogInforSunhh::change_procN( $gg{'pm'}, "$gg{'wrk_dir'}/$gg{'nprocF'}", $gg{'MAX_PROCESSES'} ); 
	my $pid = $gg{'pm'}->start and next; 
	# S1.1 : Process fastq group one by one to get _aln_pipe1.bam; ('YT' for adapter)
	for my $h1 (@{$gg{'smFq'}{$sm}{'fq_infor'}}) {
		my %fqHash = %$h1; 
		# Step 1 : Convert fastq to raw uBam : out file "$gg{'wrk_dir'}/$fqHash{'pref'}_u.bam"; 
		$gg{'doStep'}{1} and &step1_fq2uBam( $h1, "step1_cmd.fq2uBam.$fqHash{'pref'}", $gg{'wrk_dir'} ); 
		# Step 2 : Mark illumina adapters with picard MarkIlluminaAdapters; (markAdapter_bam with XT:i:##): out file "$gg{'wrk_dir'}/$fqHash{'pref'}_mrkAdp.bam"; 
		#          In this step, I change 'XT:' to 'YT:' for bwa information, because 'XT' is used in bwa; 
		$gg{'doStep'}{2} and &step2_mrkAdp( $h1, "step2_cmd.mrkAdp.$fqHash{'pref'}", $gg{'wrk_dir'} ); 
		# Step 3 : Align _mrkAdp.bam to reference fasta with BWA-MEM and merge it with uBam as aligned_pipe1.bam: out file "$gg{'wrk_dir'}/$fqHash{'pref'}_aln_pipe1.bam"; These files have all the reads and read base, and YT is recorded as adapters; 
		$gg{'doStep'}{3} and &step3_alnBam( $h1, "step3_cmd.alnBam.$fqHash{'pref'}", $gg{'wrk_dir'} ); 
	}
	# S1.2 : Merge fastq groups into single sample file. out file "$gg{'wrk_dir'}/${smID}_bySM_mrkAdp.bam"
	$gg{'doStep'}{4} and &step4_rg2sample( $gg{'smFq'}{$sm}{'fq_infor'}, "step4_cmd.rg2SM.$sm", $gg{'wrk_dir'} ); # These are aligned and coordinate-sorted files; 
	$gg{'doStep'}{5} and &step5_dedup( $sm, "step5_cmd.dedup.$sm", $gg{'wrk_dir'} ); # out file : ${smID}_bySM_dedup.bam and index; 
	$gg{'doStep'}{6} and &step6_fixBam( $sm, "step6_cmd.fixBam.$sm", $gg{'wrk_dir'} ); # out file : ${smID}_bySM_fix.bam ; 
	# S1.3 : Process single bam for each sample to get GVCFs; out file : ${sm}.g.vcf.gz
	$gg{'doStep'}{7} and &step7_bam2gvcf( $sm, "step7_cmd.bam2gvcf.$sm", $gg{'wrk_dir'} ); 
	$gg{'pm'}->finish; 
}
$gg{'pm'}->wait_all_children; 

# Section two   : Combine all GVCFs, and store their names in @{$gg{'jnGVCF_list'}}; 
### Method with CombineGVCFs : 
###   For most of gatk3 versions, I have to use CombineGVCFs instead of the following GenomicsDBImport to aggregate the GVCF files. 
###   Since the speed of combining is too slow, I want to use intervel lists to do it in parallel. 
$gg{'gvcf_list'} = [ map { "${_}.g.vcf.gz" } ( sort { $gg{'smFq'}{$a}{'in_order'} <=> $gg{'smFq'}{$b}{'in_order'} } keys %{$gg{'smFq'}} ) ]; 
$gg{'gvcf_num'}  = scalar( @{$gg{'gvcf_list'}} ); 
$cfg{'No_combineGVCF'} //= 50; 
for (my $i=0; $i<$gg{'gvcf_num'}; $i+=$cfg{'No_combineGVCF'}) {
	$gg{'gvcf_jIdx'} ++; # Joined GVCF number; 
	my @curr_glist; 
	for (my $j=$i; $j<$gg{'gvcf_num'} and $j<$i+$cfg{'No_combineGVCF'}; $j++) {
		push(@curr_glist, $gg{'gvcf_list'}[$j]); 
	}
	my $jnGVCF ; 
	if ( $gg{'doStep'}{8} ) {
		$jnGVCF = &step8_combineGVCF_interval( \@curr_glist, "$opts{'prj_ID'}_jn$gg{'gvcf_jIdx'}", '', "step8_cmd.combineGVCF.$opts{'prj_ID'}_$gg{'gvcf_jIdx'}", $gg{'wrk_dir'} ); 
	} else {
		$jnGVCF = "$opts{'prj_ID'}_jn$gg{'gvcf_jIdx'}.g.vcf.gz"; 
	}
	push(@{$gg{'jnGVCF_list'}}, $jnGVCF); 
}

# Section three : Call variants and filter data for the first time; 
### It should be possible to do the steps one by one for each scfID : 
###   call_rawV + filter_rawSNP + filter_rawInDel + Combine_filtered_SNP_and_InDel + slct_passed_data; 
$gg{'doStep'}{9} and &step9_gvcf2var( $gg{'jnGVCF_list'}, "$opts{'prj_ID'}", '', "step9_cmd.callRawV.$opts{'prj_ID'}", $gg{'wrk_dir'} ); 


################################################################################
##########    Sub-routines 
################################################################################

### OLD 
# $cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} RevertSam  I=${input_bam}  O=$wrk_dir/out_wiYT.bam  TMP_DIR=$cfg{'dir_tmp'}  ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA ATTRIBUTE_TO_CLEAR=XS SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true"; 


sub input_good {
	$gg{'usage_txt'} = <<"H1"; 
################################################################################
perl $0 fasdfasf

  -help                 Show this help. 

  -conf_file            [filename] Required. Such as 'pipe_gatk_conf'; 
                          This file tells the path information of softwares.
  -in_pref_list         [filename] Required. Not used yet. 
                          Since the input read depth determines the time cost for Haplotype Caller of GATK, it's better to sort input samples from the deepest to the least. 
                          Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2 [ <tab> PL <tab> PU <tab> Others]
                                   SM                RG_ID                 RG_LB              ReadOutPref      filename     filename  [    PL       PU       text]
                                   Repeat            Repeat                Repeat             Unique           Unique       Unique    [    Rep      Rep      Rep]
                                   GS109             GS109_Time1           GS109_Time1        GS109_Time1      W1_R1.fq.gz  W1_R2.fq.gz    illumina H0V2RADXX.1 NA
                          Example of Others : 
                                   DS=GS109_WDM_Time1;WDM_L1_I121.R[12].clean.fastq.gz;Len_100_100;;DT=2014-02-01;;CN=BFC2013288
                          'NA' means no input assigned. 
  -prj_ID               [String] Output prefix (default 'outGATK') for the whole project to merge GVCF files. Cannot contain directory; 

  -doStep               [String] Default is 'cmd', and could be 'cmd/all/1/2/3/1,2,3/3-5'; 
                                 'cmd' means only output command lines; 
                                 'all' means run all processes one by one; 
                                 Others means run the steps designed; 


  -wrk_dir              [dirname] Assign temporary directory to work in. 
  -cpuN                 [1]       Number of samples to be processed in parallel. 
  -plCatVar             [Boolean] If given, I will use perl instead of 'CombineVariants' to join the interval-called GVCFs. 
                                  The 'CatVariants' is not available in my GATK version. 

  ### Specific parameters : These parameters will over-write conf_file. 
  For HaplotypeCaller : 
    -ERC                [string] Could be 'BP_RESOLUTION' or 'GVCF'; Default BP_RESOLUTION. 
  For CombineGVCFs
    -intervalLen        [number] Default -1. If this is bigger than 0, I will combine GVCFs with interval list with multi-threads. 
  For GenotypeGVCFs     
    -CallByScf          [Boolean] If given, I'll try to call variants for each scaffold with multi-threads; 

################################################################################
H1
	defined $opts{'help'} and &LogInforSunhh::usage($gg{'usage_txt'}); 

	defined $opts{'conf_file'} or return 0; 
	defined $opts{'in_pref_list'} or return 0; 
	$opts{'prj_ID'} //= 'outGATK'; 
	$opts{'doStep'} //= 'cmd'; 
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
	### For -ERC in HaplotypeCaller : BP_RESOLUTION costs 2-3 times disk size of result file than GVCF, and the time costs are the similar. 
	###  The disk space usage depends on the sequence-depth and nucleotide divergence between the sample and the reference. 
	###  The disk space cost of ERC_BP_RES is approximately fixed for each sample. 
	###  By using BP_RES gvcf, the combined gvcf will also be BP_RES. 
	$gg{'para'}{'ERC'} = 'BP_RESOLUTION'; 
	defined $cfg{'ERC'} and $gg{'para'}{'ERC'} = $cfg{'ERC'}; 
	defined $opts{'ERC'} and $gg{'para'}{'ERC'} = $opts{'ERC'}; 
	$gg{'para'}{'ERC'} =~ m!^(GVCF|BP_RESOLUTION)$! or &stopErr("[Err] ERC ($gg{'para'}{'ERC'}) must be one of !^(GVCF|BP_RESOLUTION)\$!\n"); 
	# Get the path for diction of ref_fasta; 
	$cfg{'ref_fasta'} = &fileSunhh::_abs_path_4link( $cfg{'ref_fasta'} ); # Use the full path of 'ref_fasta'; 
	&generate_ref_idx( $cfg{'ref_fasta'} ) or &stopErr("[Err] Cannot generate enough index files for ref_fasta: $cfg{'ref_fasta'}\n"); 

	# Set interval length : 
	$gg{'para'}{'intervalLen'} = -1; 
	defined $cfg{'intervalLen'} and $gg{'para'}{'intervalLen'} = $cfg{'intervalLen'}; 
	defined $opts{'intervalLen'} and $gg{'para'}{'intervalLen'} = $opts{'intervalLen'}; 
	$gg{'para'}{'CallByScf'} = 0; 
	defined $cfg{'CallByScf'} and $gg{'para'}{'CallByScf'} = $cfg{'CallByScf'}; 
	$opts{'CallByScf'} and $gg{'para'}{'CallByScf'} = 1; 
	$gg{'para'}{'CallByScf'} =~ s!^\s*(false|F)\s*$!0!i; # True/False/0/1; 
	$gg{'para'}{'CallByScf'} =~ s!^\s*(true|T)\s*$!1!i; # True/False/0/1; 

	# Set steps to be run; 
	$gg{'stepNum'} = [1 .. 100]; 

	# Get fastq related information from in_pref_list; 
	$gg{'fq_infor'} = &load_prefList( $opts{'in_pref_list'} ); # Read in fastq information for bam output. 
	$gg{'smFq'}     = &groupFqBySM( $gg{'fq_infor'} ); 
	for my $t (split(/,/, $opts{'doStep'})) {
		$gg{'doStep'}{$t} = 1; 
	}
	defined $gg{'doStep'}{'all'} and do { map { $gg{'doStep'}{$_} = 1; } @{$gg{'stepNum'}} }; 

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

sub step1_fq2uBam {
	# Change input fastq to raw uBam files; 
	# Result files : $fqHash{'pref'}
	my ($h1, $fn_step1Cmd, $wrk_dir) = @_; 
	my %fqHash = %$h1; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	# &tsmsg("[Msg] Go to wrk_dir [$wrk_dir] from ori_dir [$ori_dir]\n"); 
	my $cmd = '' ; 

	$cmd .= "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} FastqToSam \\\n"; 
	if ( $fqHash{'fq2'} eq '' ) {
		$cmd .= "  FASTQ=$fqHash{'fq1'} \\\n"; 
	} else {
		$cmd .= "  FASTQ=$fqHash{'fq1'} FASTQ2=$fqHash{'fq2'} \\\n"; 
	}
	$cmd .= "  OUTPUT=$fqHash{'pref'}_u.bam \\\n"; 
	$cmd .= "  READ_GROUP_NAME=$fqHash{'RG'} \\\n"; 
	$cmd .= "  LIBRARY_NAME=$fqHash{'LB'} \\\n"; 
	$cmd .= "  SAMPLE_NAME=$fqHash{'SM'} \\\n"; 
	$cmd .= "  PLATFORM=$fqHash{'PL'} \\\n"; 
	$cmd .= "  SORT_ORDER=queryname \\\n"; 
	$fqHash{'CN'} ne '' and $cmd .= "  SEQUENCING_CENTER=\"$fqHash{'CN'}\" \\\n"; 
	$fqHash{'DS'} ne '' and $cmd .= "  DESCRIPTION=\"$fqHash{'DS'}\" \\\n"; 
	$fqHash{'DT'} ne '' and $cmd .= "  RUN_DATE=$fqHash{'DT'} \\\n"; 
	$fqHash{'PI'} ne '' and $cmd .= "  PREDICTED_INSERT_SIZE=$fqHash{'PI'} \\\n"; 
	$fqHash{'PG'} ne '' and $cmd .= "  PROGRAM_GROUP=$fqHash{'PG'} \\\n"; 
	$fqHash{'PM'} ne '' and $cmd .= "  PLATFORM_MODEL=$fqHash{'PM'} \\\n"; 
	$cmd .= "1> s1.std.$fqHash{'pref'} 2> s1.err.$fqHash{'pref'}\n"; 

	&fileSunhh::write2file($fn_step1Cmd, $cmd, '>'); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_step1Cmd, "# Step 1 commands: Convert fastq to raw uBam files\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_step1Cmd"); 
	chdir($ori_dir); 
	return; 
}# step1_fq2uBam() 

sub step2_mrkAdp {
	my ($h1, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my %fqHash = %$h1; 
	my $cmd = ''; 
	$cmd  = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} MarkIlluminaAdapters \\\n"; 
	$cmd .= "  I=$fqHash{'pref'}_u.bam \\\n"; 
	$cmd .= "  O=$fqHash{'pref'}_mrkAdp_wiXT.bam \\\n"; 
	$cmd .= "  M=$fqHash{'pref'}_mrkAdp_metrics.txt \\\n"; 
	$cmd .= "  TMP_DIR=$cfg{'dir_tmp'} \\\n"; 
	$fqHash{'fq2'} eq '' and $cmd .= "  MIN_MATCH_BASES_SE=6 MAX_ERROR_RATE_SE=0.2"; 
	$cmd .= "  ADAPTER_TRUNCATION_LENGTH=20 \\\n"; 
	$cmd .= "1> s2.std.$fqHash{'pref'} 2> s2.err.$fqHash{'pref'}\n"; 

	$cmd .= "$cfg{'exe_samtools'} view -h $fqHash{'pref'}_mrkAdp_wiXT.bam | \\\n"; 
	$cmd .= "  $cfg{'exe_perl'} $cfg{'pl_getSam'} -cnvt_XTi \\\n"; 
	$cmd .= "  -rawXTitag XT -tag4XTi YT | \\\n"; 
	$cmd .= "  $cfg{'exe_samtools'} view -bS -o $fqHash{'pref'}_mrkAdp.bam -\n"; 

	$cmd .= "rm -f $fqHash{'pref'}_mrkAdp_wiXT.bam\n"; 

	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 2 commands: Mark illumina adapters\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return; 
}# step2_mrkAdp() 

sub step3_alnBam {
	my ($h1, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my %fqHash = %$h1; 
	my $cmd = ''; 

	# "set -o pipefail\n\n"
	my $addDP = ''; 
	$fqHash{'fq2'} eq '' or $addDP = ( $cfg{'para_bwa'} =~ m!(^|\s)\-p(\s|$)! ) ? '' : '-p' ; 
	$cmd .= "( $cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} SamToFastq \\\n"; 
	$cmd .= "  I=$fqHash{'pref'}_mrkAdp.bam \\\n"; 
	$cmd .= "  FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=YT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true INCLUDE_NON_PRIMARY_ALIGNMENTS=false \\\n"; # Here I use 'YT' instead of 'XT' for adapter tag; 
	$cmd .= "  TMP_DIR=$cfg{'dir_tmp'} | \\\n"; 
	$cmd .= "  $cfg{'exe_bwa'} mem $cfg{'para_bwa'} $addDP \\\n"; 
	$cmd .= "  $cfg{'ref_fasta'} \\\n"; 
	$cmd .= "  /dev/stdin | \\\n"; 
	$cmd .= "  $cfg{'exe_java'} -Xmx16G -jar $cfg{'jar_picard'} MergeBamAlignment \\\n"; 
	$cmd .= "  SORT_ORDER=coordinate \\\n"; 
	# $cmd .= "  PROGRAM_RECORD_ID='bwamem' PROGRAM_GROUP_VERSION=\"bwa_version\" PROGRAM_GROUP_COMMAND_LINE=\"bwa_commandline\" PROGRAM_GROUP_NAME=\"bwamem\" \\\n"; 
	$cmd .= "  UNMAPPED_READ_STRATEGY=COPY_TO_TAG ALIGNER_PROPER_PAIR_FLAGS=true UNMAP_CONTAMINANT_READS=true MIN_UNCLIPPED_BASES=32 \\\n"; 
	$cmd .= "  IS_BISULFITE_SEQUENCE=false ALIGNED_READS_ONLY=false \\\n"; 
	if ( $fqHash{'fq2'} eq '' ) {
		$cmd .= "  CREATE_INDEX=true ADD_MATE_CIGAR=false CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \\\n"; 
	} else {
		$cmd .= "  CREATE_INDEX=true ADD_MATE_CIGAR=true  CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \\\n"; 
	}
	$cmd .= "  INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant \\\n"; 
	$cmd .= "  ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA ATTRIBUTES_TO_RETAIN=XT ATTRIBUTES_TO_RETAIN=X0 \\\n"; 
	$cmd .= "  ATTRIBUTES_TO_REMOVE=NM ATTRIBUTES_TO_REMOVE=MD \\\n"; 
	$cmd .= "  ATTRIBUTES_TO_RETAIN=YT \\\n"; 
	$cmd .= "  TMP_DIR=$cfg{'dir_tmp'} \\\n"; 
	$cmd .= "  ALIGNED_BAM=/dev/stdin   UNMAPPED_BAM=$fqHash{'pref'}_mrkAdp.bam \\\n"; 
	# $cmd .= "  ALIGNED_BAM=/dev/stdin   UNMAPPED_BAM=$fqHash{'pref'}_u.bam \\\n"; 
	$cmd .= "  OUTPUT=$fqHash{'pref'}_aln_pipe1.bam \\\n"; 
	$cmd .= "  R=$cfg{'ref_fasta'} \\\n"; 
	$cmd .= "  ) \\\n"; 
	$cmd .= "1> s3.std.$fqHash{'pref'} 2> s3.err.$fqHash{'pref'}\n"; 

	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 3 commands: Align uBams and merge alignment with adapter information.\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return; 
}# step3_alnBam() 

sub step4_rg2sample {
	my ($smAR, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my @alnPipe1_bamList = map { "$_->{'pref'}_aln_pipe1.bam" } @{$smAR}; 
	my $smID = $smAR->[0]{'SM'}; 
	my $cmd = ''; 

	if ( @alnPipe1_bamList > 1 ) {
		# Input : $gg{'wrk_dir'}/$fqHash{'pref'}_aln_pipe1.bam 
		$cmd .= "$cfg{'exe_java'} -jar $cfg{'jar_picard'} MergeSamFiles \\\n"; 
		$cmd .= "  OUTPUT=${smID}_bySM_mrkAdp.bam \\\n"; 
		$cmd .= "  SORT_ORDER=coordinate \\\n"; 
		$cmd .= "  ASSUME_SORTED=true \\\n"; 
		for ( @alnPipe1_bamList ) {
			$cmd .= "  INPUT=$_ COMMENT=Add_$_ \\\n"; 
		}
	} else {
		$cmd .= "cp -p $alnPipe1_bamList[0] ${smID}_bySM_mrkAdp.bam \\\n"; 
	}
	$cmd .= "1> s4.std.${smID} 2> s4.err.${smID}\n"; 

	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 4 commands: Merge RGs into single sample file.\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return; 
}# step4_rg2sample() 

sub step5_dedup {
	# out file : ${smID}_bySM_dedup.bam and index; 
	my ($smID, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my $cmd = ''; 

	$cmd .= "$cfg{'exe_java'} -Xmx32G -jar $cfg{'jar_picard'} MarkDuplicates \\\n"; 
	$cmd .= "  INPUT=${smID}_bySM_mrkAdp.bam \\\n"; 
	$cmd .= "  OUTPUT=${smID}_bySM_dedup.bam \\\n"; 
	$cmd .= "  METRICS_FILE=${smID}_bySM_dedup_metrics.txt \\\n"; 
	$cmd .= "  OPTICAL_DUPLICATE_PIXEL_DISTANCE=250 CREATE_INDEX=true \\\n"; # 2500 for HiSeq X, 250 for others; 
	$cmd .= "  ASSUME_SORT_ORDER=coordinate \\\n"; 
	$cmd .= "  TMP_DIR=$cfg{'dir_tmp'} \\\n"; 
	$cmd .= "1> s5.std.$smID 2> s5.err.$smID\n"; 
	
	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 5 commands: Mark duplications in sample file.\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return; 
}# step5_dedup() 

sub step6_fixBam {
	# To do : BQSR , realignment, fix SetNmMdAndUqTags; 
	my ($smID, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my $cmd = ''; 

	# realign : https://software.broadinstitute.org/gatk/documentation/article?id=7156#section1 
	# $cmd .= "$cfg{'exe_java'} -jar $cfg{'jar_gatk'} -T RealignerTargetCreator \\\n"; 
	# $cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
	# $cmd .= "  "; 

	$cmd .= "$cfg{'exe_java'} -Xmx16G -jar $cfg{'jar_picard'} SetNmMdAndUqTags \\\n"; # Parameter 'R=' is also required. 
	$cmd .= "  I=${smID}_bySM_dedup.bam \\\n"; 
	$cmd .= "  O=${smID}_bySM_fix.bam \\\n"; 
	$cmd .= "  R=$cfg{ref_fasta} \\\n"; 
	$cmd .= "  CREATE_INDEX=true \\\n"; 
	$cmd .= "1> s6.std.fixBam.$smID 2> s6.err.fixBam.$smID\n"; 

	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 6 commands: Final fix to prepare ready-to-use bam file.\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return; 
}# step6_fixBam() 

sub step7_bam2gvcf {
	# Result file : ${smID}.g.vcf.gz ; 
	# Need to check https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/master/PairedEndSingleSampleWf.gatk4.0.wdl
	#  for task HaplotypeCaller 
	# Ref : https://github.com/gatk-workflows/gatk3-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk3.wdl 
	# It takes about 2.5-3 hours to process a bam to GVCF, so I don't bother to group CHRs for speed. 
	my ($smID, $fn_cmd, $wrk_dir) = @_; 
	my $ori_dir = &fileSunhh::_abs_path("./"); 
	chdir($wrk_dir); 
	my $cmd = ''; 

	$cmd .= "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_gatk'} -T HaplotypeCaller \\\n"; 
	$cmd .= "  -R $cfg{'ref_fasta'} \\\n"; 
	$cmd .= "  -I ${smID}_bySM_fix.bam \\\n"; 
	$cmd .= "  --genotyping_mode DISCOVERY \\\n"; 
	# $cmd .= "  -stand_call_conf 30 -ERC $gg{'para'}{'ERC'} \\\n"; 
	$cmd .= "  -ERC $gg{'para'}{'ERC'} \\\n"; # I omit ' -stand_call_conf 30 ' for GVCF calling. # https://gatkforums.broadinstitute.org/gatk/discussion/7899/reference-implementation-pairedendsinglesamplewf-pipeline#MergeBamAlignment
	$cmd .= "  --max_alternate_alleles 3 --read_filter OverclippedRead \\\n"; 
	$cmd .= "  -nct $cfg{'nct_hapCaller'} \\\n"; 
	$cmd .= "  -o ${smID}.g.vcf.gz \\\n"; 
	$cmd .= "1> s7.std.bam2gvcf.$smID 2> s7.err.bam2gvcf.$smID\n"; 

	&fileSunhh::write2file( $fn_cmd, $cmd, '>' ); 
	if ($gg{'cmd'}) {
		&stdout_file( $fn_cmd, "# Step 7 commands: Call GVCFs with HaplotypeCaller from single-sample bam file.\n" ); 
		chdir($ori_dir); 
		return; 
	}
	&exeCmd_1cmd("perl $cfg{'pl_batchRun'} $fn_cmd"); 
	chdir($ori_dir); 
	return; 
}# step7_bam2gvcf() 

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
		# @seqLen = sort { $b->[1] <=> $a->[1] } @seqLen; # I can't sort sequences because the order matters in picard and GATK. 
		my @subPref_list; 
		my @toRun_intervals; # ( [[[id,start,end],[],...], interval_total_length], [], ... )
		for (my $i=0; $i<@seqLen; $i++) {
			my ($id, $len) = @{$seqLen[$i]}; 
			for (my $sP=1; $sP<$len; $sP+=$gg{'para'}{'intervalLen'}) {
				my $eP = $sP + $gg{'para'}{'intervalLen'} - 1; 
				$eP > $len and $eP = $len; 
				if (@toRun_intervals == 0) {
					push(@toRun_intervals, [ [[$id,$sP, $eP]], $eP-$sP+1 ]); 
				} elsif ( $toRun_intervals[-1][1] >= $gg{'para'}{'intervalLen'} ) {
					push(@toRun_intervals, [ [[$id,$sP, $eP]], $eP-$sP+1 ]); 
				} else {
					push(@{$toRun_intervals[-1][0]}, [$id, $sP, $eP]); 
					$toRun_intervals[-1][1] += ($eP-$sP+1); 
				}
			}
		}
		&tsmsg("[Rec] Total " . scalar(@toRun_intervals) . " intervals will be processed parallelely in Step 8 combineGVCF.\n"); 
		for (my $i=0; $i<@toRun_intervals; $i++) {
			my $sub_gvcf_pref = "${opref}.s8.$i.interval_result"; 
			push(@subPref_list, $sub_gvcf_pref); 
			my $ti_fn = "${opref}.s8.$i.interval_list"; 
			&fileSunhh::write2file($ti_fn, "$header_txt", '>'); 
			for my $iL (@{$toRun_intervals[$i][0]}) {
				&fileSunhh::write2file($ti_fn, join("\t", $iL->[0], $iL->[1], $iL->[2], '+', join("_", $iL->[0], $iL->[1], $iL->[2]))."\n", '>>'); # It doesn't matter if strand is '+'/'-' for CombineGVCF; 
			}
			$gg{'MAX_PROCESSES'} = &LogInforSunhh::change_procN( $gg{'pm'}, "$gg{'nprocF'}", $gg{'MAX_PROCESSES'} ); 
			my $pid = $gg{'pm'}->start and next; 
			&step8_combineGVCF_interval( $glist, $sub_gvcf_pref, $ti_fn, "${fn_cmd}.interval_cmd.$i", "./" ); 
			$gg{'pm'}->finish; 
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
			# $cmd .= "  -nt 20 \\\n"; 
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
		# @seqLen = sort { $b->[1] <=> $a->[1] } @seqLen; 
		my @subPref_list; 
		my @toRun_intervals; # ( [[[id, len],[],..], interval_total_length] , [], ... )
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

