#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"inBam:s", # 
	"outCnt:s", 
	"gff:s", 
	"para_htseqC:s", 
	"exe_htseqC:s", 
	"help!", 
); 

$opts{'exe_htseqC'} //= 'htseq-count'; 
$opts{'para_htseqC'} //= ' -s reverse -r name -m union --nonunique none -f bam --type CDS --idattr Parent '; 
# https://gif.biotech.iastate.edu/rnaseq-analysis-walk-through
# HTSeq2 : http://htseq.readthedocs.io/en/release_0.10.0/
# CMD: htseq-count -r pos -f bam RNA-Seq/aligned/SRR3589957_sorted.bam reference/gencode.v26lift37.annotation.sorted.gtf > SRR3589957.count
# CMD: htseq-count -s no -r pos -f bam RNA-Seq/aligned/SRR35899${i}_sorted.bam reference/gencode.v26lift37.annotation.sorted.gtf > RNA-Seq/matrix/SRR35899${i}.count 2> RNA-Seq/matrix/SRR35899${i}.log
# CMD: htseq-count -s reverse -r pos -m union --nonunique none -f bam --type CDS --idattr Parent input_srt.bam ../db/wm97pbV2ID_evm_feiID_clean.byChr.gff3 1> input_srt.bam.cnt 2> input_srt.bam.cntLog


my $help_txt = <<HH; 
################################################################################
# perl $0   -inBam hisat2_fixNH.bam -outCnt hisat2_fixNH.cnt  -gff ref_gene.gff3
#
#   Count reads in genes with htseq-count, and keep log-info in hisat2_fixNH.cnt.Log
#
# -help 
#
# -exe_htseqC     ['$opts{'exe_htseqC'}']
# -para_htseqC    ['$opts{'para_htseqC'}']
################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'inBam'}  or &LogInforSunhh::usage($help_txt); 
defined $opts{'outCnt'} or &LogInforSunhh::usage($help_txt); 
defined $opts{'gff'} or &LogInforSunhh::usage($help_txt); 

my $cmd = ""; 
$cmd .= "$opts{'exe_htseqC'} "; 
$cmd .= " $opts{'para_htseqC'} "; 
$cmd .= " $opts{'inBam'} "; 
$cmd .= " $opts{'gff'} "; 
$cmd .= " 1> $opts{'outCnt'} 2> $opts{'outCnt'}.Log"; 
&exeCmd_1cmd($cmd); 


