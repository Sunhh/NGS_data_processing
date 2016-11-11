#!/usr/bin/perl -w
# 2016-02-12 : Pipe GATK step by step for multiple input fastq files; 
# data processing steps for DNA sequencing I want to use : 
### Step 1. fastq to uBam 
### Step 2. mark illumina adapters with picard MarkIlluminaAdapters (markAdapter_bam with XT:i:##); 
### Step 3. align markAdapter_bam bam to reference fasta with BWA-MEM and merge it with uBAM as aligned_pipe1_bam; 
### Step 4. Mark duplicates in aligned_pipe1_bam as dedup_bam ; 
###### No need to realign indels because we are going to use HaplotypeCaller (assembly-based); 
###### No BQSR step included because it needs additional data and steps; (Maybe I can make one for pumpkin when having time)

### Step 5. Merge markAdapter_bam into dedup_bam as dedup_pipe2_bam (With 'YT:i:##' tag for adapter free regions)
### Step 5b. Merge bam files according to SM name. 
###### Please note that here the 'YT:i:##' corresponds to the original reads orientation, not the one in alignments. 
###### This dedup_pipe2_bam is the one I want to keep for storing data. 
### Step 6. Call GVCF files with each dedup_pipe2_bam using 'HaplotypeCaller' ; (prefix.g.vcf)
### Step 7. Combine GVCF files with each 200 single GVCF files. (*_jnX.g.vcf)
### Step 8. Call variants in combined GVCFs. (rawV.vcf)
### Step 9. Select SNP variants and do hard filtering; 
### Step 10. Select indel variants and do hard filtering; 
### Step 11. Do VQSR with both filtered SNP and indel on rawV.vcf . 
### Step 12. .... 
# 2016-02-15 Add merge sample step (step 5b). 


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
	  "singleFq!",    # Default treating input as paired-end reads. Provide this tag if processing single-end fastq. 
	"forGBS!",        # Omit '-an DP' filter if given. 
	"prj_ID:s",       # 'outGATK'. A project ID for merging GVCF files. 
	
	"do_step1!",      # Generate commands for step 1 if given. 
	"do_step2!",      # Generate commands for step 2 if given. 
	"do_step3!",      # Generate commands for step 3 if given. 
	"do_step4!",      # Generate commands for step 4 if given. 
	"do_step5!",      # Generate commands for step 5 if given. 
	"do_step5b!",     # Generate commands for step 5b if given. 
	  "step5b_o_pref_list:s", # Output a new pref_list for the next steps. 
	"do_step5c!",     # 
	"do_step6!",      # Generate commands for step 6 if given. 
	"do_step7!",      # Generate commands for step 7 if given. 
	"do_step8!",      # Generate commands for step 8 if given. 
	"do_step9!",      # Generate commands for step 9 if given. 
	"do_step10!",     # Generate commands for step 10 if given. 
	"do_step11!",     # Generate commands for step 11 if given. 
	  "step11_sampleN:i", # If sampleN < 10, I won't use 'InbreedingCoeff'. 
	
	"conf_file:s",    # Required. This file tells the path information of softwares. 
); 

################################################################################
##########    Setting basic parameters 
################################################################################
my $usage_txt = <<HH; 

perl $0 fasdfasf

	"in_pref_list:s", # Required. 
	                    Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2
	  "singleFq!",    # Default treating input as paired-end reads. Provide this tag if processing single-end fastq. 
	"forGBS!",        # Omit '-an DP' filter if given. 
	"prj_ID:s",       # 'outGATK'. A project ID for merging GVCF files. 
	
	"do_step1!",      # Generate commands for step 1 : fastq to uBam. 
	"do_step2!",      # Generate commands for step 2 : mark illumina adapters with picard MarkIlluminaAdapters (markAdapter_bam with XT:i:##); 
	"do_step3!",      # Generate commands for step 3 : align markAdapter_bam bam to reference fasta with BWA-MEM and merge it with uBAM as aligned_pipe1_bam; 
	"do_step4!",      # Generate commands for step 4 : Mark duplicates in aligned_pipe1_bam as dedup_bam ; 
	"do_step5!",      # Generate commands for step 5 : Merge markAdapter_bam into dedup_bam as dedup_pipe2_bam (With 'YT:i:##' tag for adapter free regions)
	"do_step5b!",     # Merge bam files by SM ID. 
	  "step5b_o_pref_list:s", # Output a new pref_list for the next steps. 
	"do_step5c!",     # Do BQSR for all bam files. 
	"do_step6!",      # Generate commands for step 6 : Call GVCF files with each dedup_pipe2_bam using 'HaplotypeCaller' ; (prefix.g.vcf)
	"do_step7!",      # Generate commands for step 7 : Combine GVCF files with each 200 single GVCF files. (*_jnX.g.vcf)
	"do_step8!",      # Generate commands for step 8 : Call variants in combined GVCFs. (rawV.vcf)
	"do_step9!",      # Generate commands for step 9 : Select SNP variants and do hard filtering; 
	"do_step10!",     # Generate commands for step 10 : Select indel variants and do hard filtering; 
	"do_step11!",     # Generate commands for step 11 : Do VQSR with both filtered SNP and indel on rawV.vcf . 
	
	"conf_file:s",    # Required. This file tells the path information of softwares. 

HH
&input_good() or &LogInforSunhh::usage($usage_txt); 

$opts{'prj_ID'} //= 'outGATK'; 

my %cfg; 
$cfgs_obj->getConfig('cfg_file'=>$opts{'conf_file'}, 'replace'=>1, 'hash_r'=>\%cfg); 

$opts{'do_step1'}  and mk_cmd_for_step1(); 
$opts{'do_step2'}  and mk_cmd_for_step2(); 
$opts{'do_step3'}  and mk_cmd_for_step3(); 
$opts{'do_step4'}  and mk_cmd_for_step4(); 
$opts{'do_step5'}  and mk_cmd_for_step5(); 
$opts{'do_step5b'}  and mk_cmd_for_step5b(); 
$opts{'do_step5c'}  and mk_cmd_for_step5c(); 
$opts{'do_step6'}  and mk_cmd_for_step6(); 
$opts{'do_step7'}  and mk_cmd_for_step7(); 
$opts{'do_step8'}  and mk_cmd_for_step8(); 
$opts{'do_step9'}  and mk_cmd_for_step9(); 
$opts{'do_step10'} and mk_cmd_for_step10(); 
$opts{'do_step11'} and mk_cmd_for_step11(); 

################################################################################
##########    Examples 
################################################################################
my %example; 
$example{'file_cfg'} = <<'CFG'; 

pl_getSam             /home/Sunhh/tools/github/NGS_data_processing/Proc_Sam/get_required_sam.pl

exe_perl              /usr/bin/perl
exe_java              /usr/bin/java
jar_picard            /home/Sunhh/bin/picard.jar
jar_gatk              /home/Sunhh/bin/GenomeAnalysisTK.jar
exe_samtools          /home/Sunhh/bin/samtools_1.3
exe_bwa               /home/Sunhh/bin/bwa_0710
para_bwa              -M -t 7

ref_fasta             /data/Sunhh/temp/GBS/01.cucurbita_redo/db/P1Genom_V1.scf.fa

dir_tmp               /tmp
CFG

# S1  : Result files : $th1{'pref'}_u.bam
# S2  : Result files : $th1{'pref'}_mrkAdp.bam
# S3  : Result files : $th1{'pref'}_aln_pipe1.bam
# S4  : Result files : $th1{'pref'}_dedup_pipe1.bam
# S5  : Result files : $th1{'pref'}_dedup_pipe2.bam # This is the final bam file to be kept. 
# S5b : Result files : $th1{'SM'}_merged_pipe2.bam # 
# S6  : Result files : $th1{'pref'}.g.vcf
# S7  : Result files : $opts{'prj_ID'}_jnX.g.vcf 
# S8  : Result files : $opts{'prj_ID'}_rawV.vcf
# S9  : Result files : $opts{'prj_ID'}_rawV_snps.vcf and $opts{'prj_ID'}_filt00V_snps.vcf
# S10  : Result files : $opts{'prj_ID'}_rawV_indels.vcf and $opts{'prj_ID'}_filt00V_indels.vcf
# S11 : Result files : $opts{'prj_ID'}_recal_var.vcf


################################################################################
##########    Sub-routines 
################################################################################

# Result files : $opts{'prj_ID'}_recal_var.vcf
sub mk_cmd_for_step11 {
	my %param = @_; 
	$param{'in_rawVCF'}       //= "$opts{'prj_ID'}_rawV.vcf"; 
	$param{'in_snpFiltVCF'}   //= "$opts{'prj_ID'}_filt00V_snps.vcf"; 
	$param{'in_indelFiltVCF'} //= "$opts{'prj_ID'}_filt00V_indels.vcf"; 
	$param{'step11_sampleN'}  //= $opts{'step11_sampleN'}; 
	$param{'step11_sampleN'}  //= 'NA'; # InbreedingCoeff
	my $cmd; 
	# For snps 
	$cmd = "$cfg{'exe_java'} -Xmx4G -jar $cfg{'jar_gatk'} -T VariantRecalibrator   -R $cfg{'ref_fasta'}   -input $param{'in_rawVCF'}   -recalFile $opts{'prj_ID'}_snps.recal   -tranchesFile $opts{'prj_ID'}_snps.tranches   -nt 4   -resource:$opts{'prj_ID'}_hd_snp,known=false,training=true,truth=true,prior=10.0 $param{'in_snpFiltVCF'} "; 
	$cmd .= "  -an QD "; 
	$cmd .= "  -an MQ "; 
	$cmd .= "  -an MQRankSum "; 
	$cmd .= "  -an ReadPosRankSum "; 
	$cmd .= "  -an FS "; 
	$cmd .= "  -an SOR "; 
	$opts{'forGBS'} or $cmd .= "  -an DP "; 
	($param{'step11_sampleN'} eq 'NA' or $param{'step11_sampleN'} >= 10) and $cmd .= "  -an InbreedingCoeff "; 
	$cmd .= "  -mode SNP "; 
	$cmd .= "  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "; 
	$cmd .= "  -rscriptFile $opts{'prj_ID'}_SNP_plots.R "; 
	$cmd .= "  1>s11.std.$opts{'prj_ID'}.s1 2>s11.err.$opts{'prj_ID'}.s1"; 
	print STDOUT "$cmd\n"; 
	$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_gatk'} -T ApplyRecalibration   -R $cfg{'ref_fasta'}   -input $param{'in_rawVCF'}   -mode SNP "; 
	$cmd .= "  --ts_filter_level 99.0 "; 
	$cmd .= "  -recalFile $opts{'prj_ID'}_snps.recal "; 
	$cmd .= "  -tranchesFile $opts{'prj_ID'}_snps.tranches"; 
	$cmd .= "  -o $opts{'prj_ID'}_recal_snps_raw_indels.vcf"; 
	$cmd .= "  1>s11.std.$opts{'prj_ID'}.s2 2>s11.err.$opts{'prj_ID'}.s2"; 
	print STDOUT "$cmd\n\n"; 
	# For indels 
	$cmd = "$cfg{'exe_java'} -Xmx4G -jar $cfg{'jar_gatk'} -T VariantRecalibrator   -R $cfg{'ref_fasta'}   -input $opts{'prj_ID'}_recal_snps_raw_indels.vcf   -recalFile $opts{'prj_ID'}_indels.recal   -tranchesFile $opts{'prj_ID'}_indels.tranches   -nt 4   -resource:$opts{'prj_ID'}_hd_indel,known=false,training=true,truth=true,prior=10.0 $param{'in_indelFiltVCF'} "; 
	$cmd .= "  -an QD "; 
	$opts{'forGBS'} or $cmd .= "  -an DP "; 
	$cmd .= "  -an FS "; 
	$cmd .= "  -an SOR "; 
	$cmd .= "  -an MQRankSum "; 
	$cmd .= "  -an ReadPosRankSum "; 
	($param{'step11_sampleN'} eq 'NA' or $param{'step11_sampleN'} >= 10) and $cmd .= "  -an InbreedingCoeff "; 
	$cmd .= "  -mode INDEL "; 
	$cmd .= "  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 "; 
	$cmd .= "  --maxGaussians 4 "; 
	$cmd .= "  -rscriptFile $opts{prj_ID}_INDEL_plots.R "; 
	$cmd .= "  1>s11.std.$opts{'prj_ID'}.s3 2>s11.err.$opts{'prj_ID'}.s3"; 
	print STDOUT "$cmd\n"; 
	$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_gatk'} -T ApplyRecalibration   -R $cfg{'ref_fasta'}   -input $opts{'prj_ID'}_recal_snps_raw_indels.vcf   -mode INDEL   --ts_filter_level 99.0   -recalFile $opts{'prj_ID'}_indels.recal   -tranchesFile $opts{'prj_ID'}_indels.tranches   -o $opts{'prj_ID'}_recal_var.vcf"; 
	$cmd .= "  1>s11.std.$opts{'prj_ID'}.s4 2>s11.err.$opts{'prj_ID'}.s4"; 
	print STDOUT "$cmd\n"; 
	
}# mk_cmd_for_step11() 

# Result files : $opts{'prj_ID'}_rawV_indels.vcf and $opts{'prj_ID'}_filt00V_indels.vcf
sub mk_cmd_for_step10 {
	# Ref : http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set 
	#  Here answers to 'Jeegar' tells something about 'Hard filter' on a non-human organism for VQSR. 
	#  A good place to start is by plotting the various annotations and seeing how the distributions look.
	my $cmd; 
	$cmd = "$cfg{'exe_java'} -Xmx30G -jar $cfg{'jar_gatk'} -T SelectVariants   -R $cfg{'ref_fasta'}   -V $opts{'prj_ID'}_rawV.vcf   -selectType INDEL   -o $opts{'prj_ID'}_rawV_indels.vcf 1>s10.std.$opts{'prj_ID'}.s1 2>s10.err.$opts{'prj_ID'}.s1"; 
	print STDOUT "$cmd\n"; 
	$cmd = "$cfg{'exe_java'} -Xmx30G -jar $cfg{'jar_gatk'} -T VariantFiltration   -R $cfg{'ref_fasta'}   -V $opts{'prj_ID'}_rawV_indels.vcf   -o $opts{'prj_ID'}_filt00V_indels.vcf "; 
	$cmd .= "  --filterExpression \"QD < 2.0\"               --filterName \"QD_indel\" "; 
	$cmd .= "  --filterExpression \"FS > 200.0\"             --filterName \"FS_indel\" "; 
	$cmd .= "  --filterExpression \"ReadPosRankSum < -20.0\" --filterName \"ReadPosRankSum_indel\" "; 
	$cmd .= "1>s10.std.$opts{'prj_ID'}.s1 2>s10.err.$opts{'prj_ID'}.s1"; 
	print STDOUT "$cmd\n"; 
}# mk_cmd_for_step10() 

# Result files : $opts{'prj_ID'}_rawV_snps.vcf and $opts{'prj_ID'}_filt00V_snps.vcf
sub mk_cmd_for_step9 {
	# Ref : http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set 
	#  Here answers to 'Jeegar' tells something about 'Hard filter' on a non-human organism for VQSR. 
	#  A good place to start is by plotting the various annotations and seeing how the distributions look.
	my $cmd; 
	$cmd = "$cfg{'exe_java'} -Xmx30G -jar $cfg{'jar_gatk'} -T SelectVariants   -R $cfg{'ref_fasta'}   -V $opts{'prj_ID'}_rawV.vcf   -selectType SNP   -o $opts{'prj_ID'}_rawV_snps.vcf"; 
	print STDOUT "$cmd\n"; 
	$cmd = "$cfg{'exe_java'} -Xmx30G -jar $cfg{'jar_gatk'} -T VariantFiltration   -R $cfg{'ref_fasta'}   -V $opts{'prj_ID'}_rawV_snps.vcf   -o $opts{'prj_ID'}_filt00V_snps.vcf "; 
	$cmd .= "  --filterExpression \"QD < 2.0\"              --filterName \"QD_snp\" "; 
	$cmd .= "  --filterExpression \"FS > 60.0\"             --filterName \"FS_snp\" "; 
	$cmd .= "  --filterExpression \"MQ < 40.0\"             --filterName \"MQ_snp\" "; 
	$cmd .= "  --filterExpression \"MQRankSum < -12.5\"     --filterName \"MQRankSum_snp\" "; 
	$cmd .= "  --filterExpression \"ReadPosRankSum < -8.0\" --filterName \"ReadPosRankSum_snp\" "; 
	print STDOUT "$cmd\n"; 
}# mk_cmd_for_step9() 

# Result files : $opts{'prj_ID'}_rawV.vcf
sub mk_cmd_for_step8 {
	my ($jnList_aref) = @_; 
	$jnList_aref //= []; # [prjID_jn0.g.vcf, prjID_jn1.g.vcf, ...]
	unless (@$jnList_aref > 0) {
		my $maxj = &mk_cmd_for_step7(1); 
		for (my $i=0; $i<=$maxj; $i++) {
			my $jn_fn = "$opts{'prj_ID'}_jn${i}.g.vcf"; 
			push(@$jnList_aref, $jn_fn); 
		}
	}
	my $cmd; 
	$cmd = "$cfg{'exe_java'} -Xmx50G -jar $cfg{'jar_gatk'} -T GenotypeGVCFs   -R $cfg{'ref_fasta'}   -o $opts{'prj_ID'}_rawV.vcf "; 
	for (my $i=0; $i<@$jnList_aref; $i++) {
		$cmd .= "  --variant $jnList_aref->[$i] "; 
	}
	$cmd .= "\n"; 
	print STDOUT "$cmd\n"; 
}# mk_cmd_for_step8() 

# Result files : $opts{'prj_ID'}_jnX.g.vcf 
sub mk_cmd_for_step7 {
	my ($only_cnt) = @_; 
	$only_cnt //= 0; 
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	# print STDOUT "set -o pipefail\n\n"; 
	my $no_in_batch = 200; 
	defined $cfg{'No_combineGVCF'} and $cfg{'No_combineGVCF'} > 0 and $no_in_batch = $cfg{'No_combineGVCF'}; 
	my $max_j; 
	for (my $j=0; $j*$no_in_batch < @fq_infor; $j++) {
		my $cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_gatk'} -T CombineGVCFs   -R $cfg{'ref_fasta'}   -o $opts{'prj_ID'}_jn${j}.g.vcf "; 
		for (my $i=0; $i < $no_in_batch and $i+$j*$no_in_batch < @fq_infor; $i++) {
			my %th1 = %{$fq_infor[$i+$j*$no_in_batch]}; 
			$cmd .= "  -V $th1{'pref'}.g.vcf "; 
		}
		$cmd .= "1>s7.std.$opts{'prj_ID'}_jn${j} 2>s7.err.$opts{'prj_ID'}_jn${j}"; 
		$only_cnt == 0 and print STDOUT "$cmd\n"; 
		$max_j //= $j; 
		$max_j < $j and $max_j = $j; 
	}
	return ($max_j); 
}# mk_cmd_for_step7() 

# Result files : $th1{'pref'}.g.vcf
sub mk_cmd_for_step6 {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	# print STDOUT "set -o pipefail\n\n"; 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd; 
		# Maybe consider to use Queue for parallelizing. Or simply remove '-nct 5' if it causes error. 
		my $dupFilter = ''; 
		$opts{'forGBS'} and $dupFilter = '-drf DuplicateRead'; 
		$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_gatk'} -T HaplotypeCaller   $dupFilter   -R $cfg{'ref_fasta'}   -I $th1{'pref'}_dedup_pipe2.bam   --genotyping_mode DISCOVERY   -stand_emit_conf 10   -stand_call_conf 30   -ERC GVCF   -nct 5   -o $th1{'pref'}.g.vcf 1>s6.std.$th1{'pref'} 2>s6.err.$th1{'pref'}"; 
		print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step6() 

sub mk_cmd_for_step5c {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd = "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T BaseRecalibrator   -R $cfg{'ref_fasta'}   -I $th1{'pref'}_dedup_pipe2.bam $cfg{'known_dbsnp'} $cfg{'para_bqsr'} -o $th1{'pref'}_recal_data.table > s5c.$th1{'pref'}.std 2> s5c.$th1{'pref'}.err"; 
		print STDOUT "$cmd\n"; 
		$cmd = "### To run separately: $cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T BaseRecalibrator   -R $cfg{'ref_fasta'}   -I $th1{'pref'}_dedup_pipe2.bam $cfg{'known_dbsnp'} $cfg{'para_bqsr'} -BQSR $th1{'pref'}_recal_data.table   -o $th1{'pref'}.post_recal_data.table 1>s5c.$th1{'pref'}.std.s2a 2> s5c.$th1{'pref'}.err.s2a"; 
		print STDOUT "$cmd\n"; 
		$cmd = "### To run separately: $cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T AnalyzeCovariates   -R $cfg{'ref_fasta'}   -before $th1{'pref'}_recal_data.table -after $th1{'pref'}.post_recal_data.table -plots $th1{'pref'}.recal_plots.pdf 1> s5c.$th1{'pref'}.std.s2b 2> s5c.$th1{'pref'}.err.s2b"; 
		print STDOUT "$cmd\n"; 
		$cmd = "$cfg{'exe_java'} -Xmx10G -jar $cfg{'jar_gatk'} -T PrintReads   -R $cfg{'ref_fasta'}   -I $th1{'pref'}_dedup_pipe2.bam   -BQSR $th1{'pref'}_recal_data.table   -o $th1{'pref'}.recal_reads.bam 1> s5c.$th1{'pref'}.std.s3 2> s5c.$th1{'pref'}.err.s3"; 
		print STDOUT "$cmd\n"; 
		print STDOUT "\n"; 
	}
}# mk_cmd_for_step5c ()
 
# Result files : Result files : $th1{'SM'}_merged_pipe2.bam # Sometimes we need this step. 
sub mk_cmd_for_step5b {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	$opts{'step5b_o_pref_list'} //= "step5b_out_pref_list"; 
	unlink( $opts{'step5b_o_pref_list'} ); 
	# print STDOUT "set -o pipefail\n\n"; 
	my %sm_to_bam; 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		push(@{$sm_to_bam{$th1{'SM'}}}, "$th1{'pref'}_dedup_pipe2.bam"); 
	}
	my %used; 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $smID = $th1{'SM'}; 
		defined $used{$smID} and next; 
		$used{$smID} = 1; 
		
		my @bams = @{$sm_to_bam{$smID}}; 
		@bams == 0 and &stopErr("[Err] Failed to find bam files for [$smID]\n"); 
		my $cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} MergeSamFiles   OUTPUT=${smID}_merged_dedup_pipe2.bam   SORT_ORDER=coordinate   ASSUME_SORTED=true "; 
		for my $tf (@bams) {
			$cmd .= "  INPUT=$tf COMMENT=Add_$tf "; 
		}
		print STDOUT "$cmd 1>s5b.std.${smID}.s1 2>s5b.err.${smID}.s1\n"; 
		my $fq1_fn = 'fq1'; 
		my $fq2_fn = 'fq2'; 
		( defined $th1{'fq2'} and $th1{'fq2'} ne '' ) or $fq2_fn = 'NA'; 
		&fileSunhh::write2file($opts{'step5b_o_pref_list'}, join("\t", "${smID}", 'NA', 'NA', "${smID}_merged", $fq1_fn, $fq2_fn)."\n", '>>'); 
		$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} BuildBamIndex   I=${smID}_merged_dedup_pipe2.bam 1>s5b.std.${smID}.s2 2>s5b.err.${smID}.s2"; 
		print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step5b() 

# Result files : $th1{'pref'}_dedup_pipe2.bam # This is the final bam file to be kept. 
sub mk_cmd_for_step5 {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	print STDOUT "set -o pipefail\n\n"; 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd; 
		$cmd = "( $cfg{'exe_samtools'} view -h $th1{'pref'}_dedup_pipe1.bam |   $cfg{'exe_perl'} $cfg{'pl_getSam'} -add_XTi   -bam_wiXTi $th1{'pref'}_mrkAdp.bam   -tag4XTi YT |   $cfg{'exe_samtools'} view -bS -o $th1{'pref'}_dedup_pipe2.bam - ) 1>s5.std.$th1{'pref'}.s1 2>s5.err.$th1{'pref'}.s1"; 
		print STDOUT "$cmd\n"; 
		$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} BuildBamIndex   INPUT=$th1{'pref'}_dedup_pipe2.bam 1>s5.std.$th1{'pref'}.s2 2>s5.err.$th1{'pref'}.s2"; 
		print STDOUT "$cmd\n"; 
		# The following command is not needed. 
		### $cmd = "$cfg{'exe_java'} -Xmx16G -jar $cfg{'jar_picard'} MergeBamAlignment    R=$cfg{'ref_fasta'}   ALIGNED_BAM=$th1{'pref'}_dedup_pipe2.bam   UNMAPPED_BAM=$th1{'pref'}_u.bam   CREATE_INDEX=true ADD_MATE_CIGAR=true   CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true   INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1   PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS   ATTRIBUTES_TO_RETAIN=XA   ATTRIBUTES_TO_RETAIN=XT   ATTRIBUTES_TO_RETAIN=YT   TMP_DIR=$cfg{'dir_tmp'} 1>s5.std.$th1{'pref'} 2>s5.err.$th1{'pref'}"; 
		### print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step5() 

# Result files : $th1{'pref'}_dedup_pipe1.bam
sub mk_cmd_for_step4 {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	print STDOUT "set -o pipefail\n\n"; 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd; 
		$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} SortSam   INPUT=$th1{'pref'}_aln_pipe1.bam   OUTPUT=$th1{'pref'}_srt_pipe1.bam   SORT_ORDER=coordinate"; 
		print STDOUT "$cmd\n"; 
		### I don't want to use the tool 'MarkDuplicatesWithMateCigar' currently, because I don't know what problems it might cause. 
		$cmd = "$cfg{'exe_java'} -Xmx32G -jar $cfg{'jar_picard'} MarkDuplicates   INPUT=$th1{'pref'}_srt_pipe1.bam   OUTPUT=$th1{'pref'}_dedup_pipe1.bam   METRICS_FILE=$th1{'pref'}_dedup_pipe1_metrics.txt   OPTICAL_DUPLICATE_PIXEL_DISTANCE=250   CREATE_INDEX=true   TMP_DIR=$cfg{'dir_tmp'}"; 
		print STDOUT "$cmd\n"; 
		$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} BuildBamIndex   INPUT=$th1{'pref'}_dedup_pipe1.bam"; 
		print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step4() 


sub generate_bwaMEM_idx {
	my $refDb = shift; 
	my $is_good = 1; 
	for (qw/amb ann bwt pac sa/) {
		my $fn = "$refDb.$_"; 
		-e $fn or do { $is_good = 0; last; }; 
	}
	if ($is_good == 0) {
		&exeCmd_1cmd("$cfg{'exe_bwa'} index $refDb 1>$refDb.dbIdx.std 2>$refDb.dbIdx.err") and return 0; 
	}
	return 1; 
}

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
	return 1; 
}# generate_ref_idx() 

# Result files : $th1{'pref'}_aln_pipe1.bam
sub mk_cmd_for_step3 {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	print STDOUT "set -o pipefail\n\n"; 
	# Check if we have enough files for ref_fasta; 
	&generate_ref_idx($cfg{'ref_fasta'}) or &stopErr("[Err] Cannot generate enough index files for ref_fasta: $cfg{'ref_fasta'}\n"); 
	&generate_bwaMEM_idx($cfg{'ref_fasta'}) or &stopErr("[Err] Failed to generate bwaMEM-db-index for [$cfg{'ref_fasta'}]\n"); 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd; 
		if ($th1{'fq2'} eq '') {
			$cmd = "( $cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} SamToFastq   I=$th1{'pref'}_mrkAdp.bam   FASTQ=/dev/stdout   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true   TMP_DIR=$cfg{'dir_tmp'} |   $cfg{'exe_bwa'} mem $cfg{'para_bwa'}   $cfg{'ref_fasta'}   /dev/stdin |   $cfg{'exe_java'} -Xmx16G -jar $cfg{'jar_picard'} MergeBamAlignment   ALIGNED_BAM=/dev/stdin   UNMAPPED_BAM=$th1{'pref'}_u.bam   OUTPUT=$th1{'pref'}_aln_pipe1.bam   R=$cfg{'ref_fasta'}   CREATE_INDEX=true ADD_MATE_CIGAR=false   CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true   INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1   PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS   ATTRIBUTES_TO_RETAIN=XA   ATTRIBUTES_TO_RETAIN=XT   TMP_DIR=$cfg{'dir_tmp'} ) 1>s3.std.$th1{'pref'} 2>s3.err.$th1{'pref'}"; 
		} else {
			my $addDP = ( $cfg{'para_bwa'} =~ m!(^|\s)\-p(\s|$)! ) ? '' : '-p' ; 
			$cmd = "( $cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} SamToFastq   I=$th1{'pref'}_mrkAdp.bam   FASTQ=/dev/stdout   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true   TMP_DIR=$cfg{'dir_tmp'} |   $cfg{'exe_bwa'} mem $cfg{'para_bwa'}   $addDP   $cfg{'ref_fasta'}   /dev/stdin |   $cfg{'exe_java'} -Xmx16G -jar $cfg{'jar_picard'} MergeBamAlignment   ALIGNED_BAM=/dev/stdin   UNMAPPED_BAM=$th1{'pref'}_u.bam   OUTPUT=$th1{'pref'}_aln_pipe1.bam   R=$cfg{'ref_fasta'}   CREATE_INDEX=true ADD_MATE_CIGAR=true   CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true   INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1   PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS   ATTRIBUTES_TO_RETAIN=XA   ATTRIBUTES_TO_RETAIN=XT   TMP_DIR=$cfg{'dir_tmp'} ) 1>s3.std.$th1{'pref'} 2>s3.err.$th1{'pref'}"; 
		}
		print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step3() 


# Result files : $th1{'pref'}_mrkAdp.bam
sub mk_cmd_for_step2 {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd; 
		if ($th1{'fq2'} eq '') {
			$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} MarkIlluminaAdapters   I=$th1{'pref'}_u.bam   O=$th1{'pref'}_mrkAdp.bam   M=$th1{'pref'}_mrkAdp_metrics.txt   TMP_DIR=$cfg{'dir_tmp'}   MIN_MATCH_BASES_SE=6   MAX_ERROR_RATE_SE=0.2   ADAPTER_TRUNCATION_LENGTH=20 1>s2.std.$th1{'pref'} 2>s2.err.$th1{'pref'}"; 
		} else {
			$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} MarkIlluminaAdapters   I=$th1{'pref'}_u.bam   O=$th1{'pref'}_mrkAdp.bam   M=$th1{'pref'}_mrkAdp_metrics.txt   TMP_DIR=$cfg{'dir_tmp'}   ADAPTER_TRUNCATION_LENGTH=20 1>s2.std.$th1{'pref'} 2>s2.err.$th1{'pref'}"; 
		}
		print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step2() 

# Result files : $th1{'pref'}_u.bam 
sub mk_cmd_for_step1 {
	my @fq_infor = @{ &load_prefList($opts{'in_pref_list'}) }; # 
	for my $h1 (@fq_infor) {
		my %th1 = %$h1; 
		my $cmd; 
		$th1{'PL'} //= 'illumina'; 
		$th1{'PU'} //= $th1{'RG'}; 
		if ($th1{'fq2'} eq '') {
			$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} FastqToSam   FASTQ=$th1{'fq1'}   OUTPUT=$th1{'pref'}_u.bam   READ_GROUP_NAME=$th1{'RG'}   LIBRARY_NAME=$th1{'LB'}   SAMPLE_NAME=$th1{'SM'}   PLATFORM=$th1{'PL'}   PLATFORM_UNIT=$th1{'PU'}   SORT_ORDER=queryname 1>s1.std.$th1{'pref'} 2>s1.err.$th1{'pref'}"; 
		} else {
			$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} FastqToSam   FASTQ=$th1{'fq1'}   FASTQ2=$th1{'fq2'}   OUTPUT=$th1{'pref'}_u.bam   READ_GROUP_NAME=$th1{'RG'}   LIBRARY_NAME=$th1{'LB'}   SAMPLE_NAME=$th1{'SM'} PLATFORM=$th1{'PL'}   PLATFORM_UNIT=$th1{'PU'}   SORT_ORDER=queryname 1>s1.std.$th1{'pref'} 2>s1.err.$th1{'pref'}"; 
		}
		print STDOUT "$cmd\n"; 
	}
}# mk_cmd_for_step1() 



# check if the input is sufficient to go on. 
sub input_good {
	defined $opts{'in_pref_list'} or return 0; 
	defined $opts{'conf_file'} or return 0; 
	return 1; 
}# input_good() 

# I would lines heading with '#' 
# Return for paired: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>fq2, 'PL'=>'illumina', 'PU'=>RG}, {}, ... ])
# Return for single: ([ {'SM'=>SM, 'RG'=>RG, 'LB'=>LB, 'pref'=>pref, 'fq1'=>fq1, 'fq2'=>'' , 'PL'=>'illumina', 'PU'=>RG}, {}, ... ])
sub load_prefList {
	my $fn = shift; 
	my $fh = &openFH($fn, '<'); 
	# Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2 <tab> PL <tab> PU
	my @back; 
	while (<$fh>) {
		m/^\s*($|#)/ and next; 
		chomp; 
		my @ta = split(/\t/, $_); 
		my ($sm, $rg, $lb, $pref, $fq1, $fq2, $pl, $pu) = @ta; 
		$pl //= 'illumina'; 
		$pu //= $rg; 
		my %th; 
		$opts{'singleFq'} and $fq2 = ''; 
		$fq2 //= ''; 
		$fq1 =~ m/^na$/i and do { &tsmsg("[Wrn] Skip bad line without fq1 [$_]\n"); next; }; 
		$fq2 =~ m/^na$/i and $fq2 = ''; 
		$th{'SM'} = $sm; 
		$th{'RG'} = $rg; 
		$th{'LB'} = $lb; 
		$th{'pref'} = $pref; 
		$th{'fq1'} = $fq1; 
		$th{'fq2'} = $fq2; 
		$th{'PL'}  = $pl; 
		$th{'PU'}  = $pu; 
		push(@back, \%th); 
	}
	close($fn); 
	return(\@back); 
}# load_prefList() 


################################################################################
##########    Examples 
################################################################################
