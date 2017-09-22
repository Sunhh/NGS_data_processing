#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SeqAlnSunhh; 
my $sas = SeqAlnSunhh->new(); 

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"db:s", # /data/Sunhh/database/db_bwa/Watermelon/WM97_v6.chr.fa
	"dbTag:s", # to97chr
	"exe_bwa:s", # bwa
	"exe_samtools:s", # samtools
	"step2do:s@", # Steps to be done, default is ('all'). Could be one of qw/all aln_r1 aln_r2 sampe bam_sort bam_index rm_sai rm_rawbam sampe2bam sampe2sam sam2bam rm_rawsam/ for PE alignment. 
	"oBamPre:s", 
	"para_aln:s", # -t $cpuN -n 0.03 -o 1 -e 2
	"para_sampe:s", 
	"para_samse:s", 
	"para_bamSort:s", # -@ $cpuN -m $mem_limit 
	"para_sam2bam:s", 
	"printCmd!", 

	"aln_type:s", # PE
	"inFq1:s", "inFq2:s", 
	"para_list:s", 
	"help!", 
); 

sub usage {
	print <<HH;
################################################################################
#  perl $0 
#
# -para_list  [''] A file recording all information needed in bwa alignment. 
#                  This will overwrite the general settings. 
#
#  -inFq1     input_R1 fastq
#  -inFq2     input_R2 fastq
#  -db        bwa database indexed
#  -dbTag     A string tagged following the .sai files. 
#
#  -aln_type        ['PE'] Paired-end(PE) alignment or Single-end(SE) alignment. 
#
#  -exe_bwa         ['bwa'] in default 
#  -exe_samtools    ['samtools'] 
#
#  -step2do         ['all'] in default. could be assigned multiple times. qw/all aln_r1 aln_r2 sampe bam_sort bam_index rm_sai rm_rawbam sampe2bam sampe2sam sam2bam rm_rawsam/
#  
#  -oBamPre         Prefix of output files. 
#  -para_aln        [''] used in 'bwa aln'. Could be '-t \$cpuN -n 0.03 -o 1 -e 2'
#  -para_sampe      [''] used in 'bwa sampe'. Could be ' -s'
#  -para_samse      [''] used in 'bwa samse'. 
#  -para_bamSort    [''] used in 'samtools sort'. Could be '-@ 20 -m 5G' if supported. 
#  -para_sam2bam    [''] used in 'samtools view' converting sam to bam file. 
#  
#  -printCmd        [FALSE] Only output commands if given. 
################################################################################
# Example of -para_list
#oBamPre	inFq1	inFq2
#Slycopersicoides0_P	Slycopersicoides0_R1.fq.gz	Slycopersicoides0_R2.fq.gz
#Slycopersicoides1_P	Slycopersicoides1_R1.fq.gz	Slycopersicoides1_R2.fq.gz
#Slycopersicoides2_P	Slycopersicoides2_R1.fq.gz	Slycopersicoides2_R2.fq.gz
#Slycopersicoides3_P	Slycopersicoides3_R1.fq.gz	Slycopersicoides3_R2.fq.gz
#Slycopersicoides10_P	Slycopersicoides10_R1.fq.gz	Slycopersicoides10_R2.fq.gz
################################################################################
HH
	exit; 
}

$opts{'help'} and &usage(); 

for (qw/db dbTag oBamPre/) {
	$opts{$_} //= 'NA'; 
}
$opts{'exe_bwa'} //= 'bwa'; 
$opts{'exe_samtools'} //= 'samtools'; 
$opts{'step2do'} //= ['all']; 
$opts{'printCmd'} //= 0; 
$opts{'aln_type'} //= 'PE'; 

my @batch; 
if (defined $opts{'inFq1'}) {
	push(@batch, { %opts }); 
}
if ( defined $opts{'para_list'} ) {
	my @header; 
	open F,'<',"$opts{'para_list'}" or die; 
	my $hh = <F>; 
	chomp($hh); 
	@header = split(/\t/, $hh); 
	while (<F>) {
		chomp; m/^(\s*$|\s*#)/ and next; 
		my @ta = split(/\t/, $_); 
		my %cur_opts; 
		for ( my $i=0; $i<@header; $i++ ) {
			if ( $header[$i] eq 'step2do' ) {
				$ta[$i] = [ map { $_ =~ s/\s//g; $_; } split(/,/, $ta[$i]) ]; 
			}
			$cur_opts{$header[$i]} = $ta[$i]; 
		}
		push(@batch, { %cur_opts }); 
		for my $tk ( keys %opts ) {
			defined $batch[-1]{$tk} or $batch[-1]{$tk} = $opts{$tk}; 
		}
	}
	close F; 
}

$#batch == -1 and &usage(); 

for (my $i=0; $i<@batch; $i++) {
	$sas->bwaAln( %{$batch[$i]} ); 
}

&tsmsg("[Rec]All done.\n"); 
