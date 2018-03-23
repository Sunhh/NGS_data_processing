#!/usr/bin/perl
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
	"log_lineN:i", 
	"help!", 
); 

$opts{'exe_samtools'} //= 'samtools'; 
$opts{'log_lineN'}    //= 1e6; 

my $help_txt = <<HH; 
######################################################################
# perl $0 -ibam input.bam > input.bam.rdNum
#
# -exe_samtools     [$opts{'exe_samtools'}]. 
# -input_sam        [Boolean]
# -log_lineN        [$opts{'log_lineN'}]
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

print join("\t", qw/InFile  Total_size  Total_Rd_num  Mean_Rd_size  Range_Rd_size  PhredCut  Time/)."\n"; 
for my $cur_bam (@{$opts{'ibam'}}) {

	my $sam_fh; 
	my ($rdN, $rdSize, $max, $min); 
	my %h; 
	if ($opts{'input_sam'}) {
		$sam_fh = &openFH($cur_bam, '<'); 
	} else {
		$sam_fh = &SeqAlnSunhh::openSam( $cur_bam, undef(), { 'wiH'=>0, 'verbose'=>1, 'exe_samtools'=>$opts{'exe_samtools'} } );
	}

	my %cnt; 
	$cnt{'log_section'} = { 'cntN_base' => 0, 'cntN_step' => $opts{'log_lineN'} }; 
	while (<$sam_fh>) {
		&fileSunhh::log_section( $. , $cnt{'log_section'} ) and &tsmsg("[Msg][PID=$$] Reading [$cur_bam] $. line.\n");
		m!^\@! and next; 
		chomp; 
		my @ta = &splitL("\t", $_); 
		defined $flag2R{$ta[1]} or die "Failed to get R_rank for flag [$ta[1]]\n"; 
		my $k = "$ta[0]\t$flag2R{$ta[1]}"; 
		defined $h{$k} and next; 
		$h{$k} = 1; 
		my $len = length($ta[9]); 
		$rdN ++; 
		$rdSize += $len; 
		$max //= $len; 
		$min //= $len; 
		$max < $len and $max = $len; 
		$min > $len and $min = $len; 
	}
	close($sam_fh); 
	$rdN //= 0; $rdSize //= 0; $max //= 0; $min //= 0; 
	my $avgLen = ($rdN > 0) ? $rdSize/$rdN : 0 ; 
	print join("\t", $cur_bam, $rdSize, $rdN, $avgLen, "${min}-${max}", 'NA', scalar(localtime()))."\n"; 
}



