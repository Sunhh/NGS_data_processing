#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"jar_picard:s", 
		"maxDropFrac:f", 
	"exe_java:s", 
	"exe_samtools:s", 
	"inBam:s", 
	"outBam:s", 
	"outFq:s", 
); 
$opts{'exe_java'}   //= '/usr/java/jre1.8.0_144/bin/java'; 
$opts{'jar_picard'} //= "/home/Sunhh/src/align/picard/v2.10.3/picard.jar"; 
$opts{'exe_samtools'} //= '/opt/align/samtools/samtools-1.5/samtools'; 

$opts{'maxDropFrac'} //= 0.001; 

my $help_txt = <<HH; 
################################################################################
# perl $0 -inBam  BF80d2_H3C7HCCXY_aln_pipe1.bam -outBam BF80d2_H3C7HCCXY_u.bam
#
# -help 
#
# -exe_java            [$opts{'exe_java'}] 
# -jar_picard          [$opts{'jar_picard'}] 
# -exe_samtools        [$opts{'exe_samtools'}] 
#
# -maxDropFrac         [$opts{'maxDropFrac'}] For MAX_DISCARD_FRACTION; 
# 
################################################################################
HH

( defined $opts{'inBam'} and defined $opts{'outBam'} ) or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my $alnBam = $opts{'inBam'}; 
my $outBam = $opts{'outBam'}; 
my $tmp_dir = &fileSunhh::new_tmp_dir('create' => 1); 
defined $tmp_dir or &stopErr("[Err] Failed to create tmp_dir for $$\n"); 

&tsmsg("[Msg] Begin to convert $alnBam to $outBam\n"); 

# I cann't ask picard RevertSam to output to STDOUT, so I have to use a temporary file. 
my $openCmd = ''; 
$openCmd .= "$opts{'exe_java'} -Xmx8G -jar $opts{'jar_picard'} RevertSam "; 
$openCmd .= "  I=$alnBam "; 
$openCmd .= "  O=${tmp_dir}/a.bam "; # OUTPUT_BY_READGROUP=true for O=dir
$openCmd .= "  SORT_ORDER=queryname SANITIZE=true MAX_DISCARD_FRACTION=$opts{'maxDropFrac'} "; # 
$openCmd .= "  REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true RESTORE_ORIGINAL_QUALITIES=true "; 
# $openCmd .= "  ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=X0 "; 
# $openCmd .= "  ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP "; 
&runCmd( $openCmd ); 

my $outCmd = ''; 
$outCmd .= "$opts{'exe_samtools'} view -o $outBam -"; 
# open F,'-|', "$openCmd" or die "Failed to run input CMD: $openCmd\n"; 
open F,'-|', "$opts{'exe_samtools'} view -h $tmp_dir/a.bam" or die "Failed to open file $tmp_dir/a.bam\n"; 
open O,'|-', "$outCmd" or die "Failed to run output CMD: $outCmd\n"; 
while (<F>) {
	chomp; 
	if (m!^\@!) {
		print O "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	my @o = @ta[0 .. 10]; 
	for (my $i=11; $i<@ta; $i++) {
		$ta[$i] =~ m!^RG:! and push(@o, $ta[$i]); 
	}
	print O join("\t", @o)."\n"; 
}
close(F); 
close(O); 


if (defined $opts{'outFq'} and $opts{'outFq'} ne '') {
	$outCmd = ''; 
	$outCmd .= "$opts{'exe_java'} -Xmx8G -jar $opts{'jar_picard'} SamToFastq "; 
	$outCmd .= " NON_PF=true "; 
	$outCmd .= " I=$outBam "; 
	$outCmd .= " F=$opts{'outFq'}_R1.fq "; 
	$outCmd .= " F2=$opts{'outFq'}_R2.fq"; 
	&runCmd( $outCmd ); 
}

&fileSunhh::_rmtree($tmp_dir); 
&tsmsg("[Msg] Finish to convert $alnBam to $outBam\n"); 

sub runCmd {
	&exeCmd_1cmd( $_[0] ) and &stopErr("[Err] Failed at CMD: $_[0]\n"); 
	return; 
}# runCmd 

