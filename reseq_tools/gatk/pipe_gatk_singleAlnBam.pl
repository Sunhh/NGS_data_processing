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
	"prj_ID:s",       # 'outGATK'. A project ID for merging GVCF files. 
	
	"conf_file:s",    # Required. This file tells the path information of softwares. 
); 

################################################################################
##########    Setting basic parameters 
################################################################################
my $usage_txt = <<HH; 

perl $0 fasdfasf

	"in_pref_list:s", # Required. 
	                    Format : SAMPLE_NAME <tab> READ_GROUP_NAME <tab> LIBRARY_NAME <tab> dataPrefix <tab> in_fq1 <tab> in_fq2
	"prj_ID:s",       # 'outGATK'. A project ID for merging GVCF files. 
	
	"conf_file:s",    # Required. This file tells the path information of softwares. 

HH
&input_good() or &LogInforSunhh::usage($usage_txt); 

# check if the input is sufficient to go on. 
sub input_good {
	defined $opts{'conf_file'} or return 0; 
	return 1; 
}# input_good() 


$opts{'prj_ID'} //= 'outGATK'; 

my %cfg; 
$cfgs_obj->getConfig('cfg_file'=>$opts{'conf_file'}, 'replace'=>1, 'hash_r'=>\%cfg); 
$cfg{'nct_hapCaller'} //= 5; 


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

nct_hapCaller         5
CFG

################################################################################
##########    Examples 
################################################################################


################################################################################
##########    Sub-routines 
################################################################################

my $input_bam = shift; 
my $opref = $opts{'prj_ID'}; 
my $wrk_dir = &fileSunhh::new_tmp_dir('create' => 1) or &stopErr("[Err] Failed to create tmp dir for [$input_bam]\n"); 
&fileSunhh::write2file("$wrk_dir/input", "Input=$input_bam\n", '>'); 
&fileSunhh::write2file("$wrk_dir/input", "OutPref=$opref\n", '>>'); 

sub save_and_exeCmd {
	my $c = shift; 
	&fileSunhh::write2file("$wrk_dir/cmd", "$c\n", '>>'); 
	&exeCmd_1cmd($c) and &stopErr("[Err] Failed to run cmd: $c\n");
	return(); 
}# save_and_exeCmd() 
my $cmd; 

if (1) {
	$input_bam =~ m!^(\S+)\.bam$! or die "Bad input bam file name [$input_bam].\n"; 
	my $input_pref = $1; 
	unless ( -e "${input_pref}.bai" or -e "${input_pref}.bam.bai" ) {
		$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} BuildBamIndex I=$input_bam"; 
		&save_and_exeCmd($cmd); 
	}
}

$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} RevertSam  I=${input_bam}  O=$wrk_dir/out_wiYT.bam  TMP_DIR=$cfg{'dir_tmp'}  ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA ATTRIBUTE_TO_CLEAR=XS SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true"; 
&save_and_exeCmd($cmd); 

$cmd = "($cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_picard'} SamToFastq  I=$wrk_dir/out_wiYT.bam  FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=YT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true  TMP_DIR=$cfg{'dir_tmp'} | $cfg{'exe_bwa'} mem $cfg{'para_bwa'} -p $cfg{'ref_fasta'} /dev/stdin | $cfg{'exe_java'} -Xmx16G -jar $cfg{'jar_picard'} MergeBamAlignment  ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=$wrk_dir/out_wiYT.bam  OUTPUT=$wrk_dir/out_aln_pipe1.bam R=$cfg{'ref_fasta'} CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA ATTRIBUTES_TO_RETAIN=YT  TMP_DIR=$cfg{'dir_tmp'} ) 1>$wrk_dir/s3.std.out 2>$wrk_dir/s3.err.out"; 
&save_and_exeCmd($cmd);

$cmd = "$cfg{'exe_java'} -jar $cfg{'jar_picard'} SortSam INPUT=$wrk_dir/out_aln_pipe1.bam OUTPUT=$wrk_dir/out_srt_pipe1.bam SORT_ORDER=coordinate"; 
&save_and_exeCmd($cmd);
$cmd = "$cfg{'exe_java'} -Xmx32G -jar $cfg{'jar_picard'} MarkDuplicates  INPUT=$wrk_dir/out_srt_pipe1.bam OUTPUT=$wrk_dir/out_dedup_pipe1.bam METRICS_FILE=$wrk_dir/out_dedup_pipe1_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=250 CREATE_INDEX=true TMP_DIR=$cfg{'dir_tmp'}"; 
&save_and_exeCmd($cmd); 

$cmd = "$cfg{'exe_java'} -Xmx8G -jar $cfg{'jar_gatk'} -T HaplotypeCaller  -R $cfg{'ref_fasta'} -I $wrk_dir/out_dedup_pipe1.bam  --genotyping_mode DISCOVERY -stand_call_conf 30 -ERC GVCF  -nct $cfg{'nct_hapCaller'} -o $wrk_dir/out_merged.g.vcf 1>$wrk_dir/s6.std.out 2>$wrk_dir/s6.err.out"; 
&save_and_exeCmd($cmd); 

$cmd = "mv $wrk_dir/out_merged.g.vcf ${opref}.g.vcf"; 
&save_and_exeCmd($cmd); 

&save_and_exeCmd("mv $wrk_dir/out_dedup_pipe1.bam ${opref}_dedup_pipe2.bam"); 
&save_and_exeCmd("mv $wrk_dir/cmd ${opref}.cmd"); 
$cmd = "gzip ${opref}.g.vcf"; 
&save_and_exeCmd($cmd); 

&fileSunhh::_rmtree($wrk_dir); 

