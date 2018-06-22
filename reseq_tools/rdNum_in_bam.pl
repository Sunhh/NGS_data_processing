#!/usr/bin/perl
# 2018-06-22 
use strict; 
use warnings; 
use SeqAlnSunhh;
use fileSunhh;
use LogInforSunhh; 
use Getopt::Long;
my %opts;
GetOptions(\%opts,
	"ibam:s@", 
	"input_sam!", 
	"exe_samtools:s",     # samtools_1.3 
	"revertBam!", "exe_java:s", "jar_picard:s", 
	"log_lineN:i", 
	"help!", 
); 

$opts{'exe_samtools'} //= 'samtools'; 
# $opts{'exe_samtools'} //= '/opt/align/samtools/samtools-1.5/samtools';
$opts{'log_lineN'}    //= 1e6; 
$opts{'exe_java'}   //= '/usr/java/jre1.8.0_144/bin/java';
$opts{'jar_picard'} //= "/home/Sunhh/src/align/picard/v2.10.3/picard.jar";


my $help_txt = <<HH; 
######################################################################
# perl $0 -ibam input.bam > input.bam.rdNum
#
# -exe_samtools     [$opts{'exe_samtools'}]. 
# -input_sam        [Boolean]
# -log_lineN        [$opts{'log_lineN'}]
#
# -revertBam        [Boolean] If this is given, I need -java_picard to convert aligned_bam to unaligned_bam at first. 
#                             This is not necessary if the aligned_bam is stored by my merged uBam pipeline. 
#   -exe_java       [$opts{'exe_java'}]
#   -jar_picard     [$opts{'jar_picard'}]
#
HH

defined $opts{'ibam'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %flag_aln_R1 = %{ &SeqAlnSunhh::mk_flag( 'keep' => '6=1,7=0' ) }; 
my %flag_aln_R2 = %{ &SeqAlnSunhh::mk_flag( 'keep' => '6=0,7=1' ) }; 
my %flag2R; 
for ( keys %flag_aln_R1 ) {
	$flag2R{$_} = 'R1'; 
}
for ( keys %flag_aln_R2 ) {
	$flag2R{$_} = 'R2'; 
}

my $tmp_dir = &fileSunhh::new_tmp_dir('create'=>1); 

print join("\t", qw/InFile  Total_size  Total_Rd_num  Mean_Rd_size  Range_Rd_size  PhredCut  Time/)."\n"; 
for my $cur_bam (@{$opts{'ibam'}}) {

	my $sam_fh; 
	if ( $opts{'revertBam'} ) {
		my $openCmd = ''; 
		$openCmd .= "$opts{'exe_java'} -Xmx8G -jar $opts{'jar_picard'} RevertSam "; 
		$openCmd .= "  I=$cur_bam "; 
		$openCmd .= "  O=${tmp_dir}/a.bam "; # OUTPUT_BY_READGROUP=true for O=dir 
		$openCmd .= "  SORT_ORDER=queryname SANITIZE=true MAX_DISCARD_FRACTION=0.001 "; #
		$openCmd .= "  REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true RESTORE_ORIGINAL_QUALITIES=true ";
		&exeCmd_1cmd( $openCmd ); 
		$sam_fh = &SeqAlnSunhh::openSam( "$tmp_dir/a.bam", undef(), { 'wiH'=>0, 'verbose'=>1, 'exe_samtools'=>$opts{'exe_samtools'} } ); 
	} elsif ( $opts{'input_sam'} ) {
		$sam_fh = &openFH($cur_bam, '<'); 
	} else {
		$sam_fh = &SeqAlnSunhh::openSam( $cur_bam, undef(), { 'wiH'=>0, 'verbose'=>1, 'exe_samtools'=>$opts{'exe_samtools'} } );
	}

	my %cnt; 
	$cnt{'log_section'} = { 'cntN_base' => 0, 'cntN_step' => $opts{'log_lineN'} }; 
	&cnt_inFh( $sam_fh, \%cnt ); 
	print join("\t", $cur_bam, $cnt{'rdSize'}, $cnt{'rdN'}, $cnt{'avgLen'}, "$cnt{'min'}-$cnt{'max'}", 'NA', scalar(localtime()))."\n"; 
}
$opts{'revertBam'} and &fileSunhh::_rmtree($tmp_dir); 

sub cnt_inFh {
	my ($sam_fh, $cntH) = @_; 
	my %h; 
	while (<$sam_fh>) {
		&fileSunhh::log_section( $. , $cntH->{'log_section'} ) and &tsmsg("[Msg][PID=$$] Reading $. line.\n"); 
		m!^\@! and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		unless ($opts{'revertBam'}) {
			defined $flag2R{$ta[1]} or die "Failed to get R_rank for flag [$ta[1]]\n"; 
			my $k = "$ta[0]\t$flag2R{$ta[1]}"; 
			defined $h{$k} and next; 
			$h{$k} = 1; 
		}
		my $len = length($ta[9]); 
		$cntH->{'rdN'} ++; 
		$cntH->{'rdSize'} += $len; 
		$cntH->{'max'} //= $len; 
		$cntH->{'min'} //= $len; 
		$cntH->{'max'} < $len and $cntH->{'max'} = $len; 
		$cntH->{'min'} > $len and $cntH->{'min'} = $len; 
	}
	close($sam_fh); 
	for (qw/rdN rdSize max min/) {
		$cntH->{$_} //= 0; 
	}
	$cntH->{'avgLen'} = ( $cntH->{'rdN'} > 0 ) ? $cntH->{'rdSize'}/$cntH->{'rdN'} : 0 ; 
	return($cntH); 
}# cnt_inFh() 


