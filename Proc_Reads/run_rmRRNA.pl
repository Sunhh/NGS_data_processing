#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 inLis\n"; 

my $dbRRNA = '/share1/db_bowtie/rRNA_silva111'; 
my $cpuN = 10; 

my $pl_extractFq = '/home/Sunhh/tools/github/NGS_data_processing/extract_fq_by_list.pl'; 


while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $pref = "$ta[0]"; 
	my $fqFile = "${pref}_highQ.single"; 
	my $oFqFile="${pref}_rmRRNA.fq"; 
	&tsmsg("[Rec] Dealing with [$fqFile]\n"); 
	&exeCmd("bowtie -v 3 -k 1 -S -p $cpuN $dbRRNA $fqFile | samtools view -S -F 4 -hb -o $fqFile.bam -"); 
	&exeCmd("samtools view $fqFile.bam | cut -f 1 > $fqFile.bam.rd"); 
	&exeCmd("perl $pl_extractFq -mode drop -rdKey -refLis $fqFile.bam.rd -srcFq $fqFile -outFq $oFqFile"); 
}
&tsmsg("[Rec] All done.\n"); 
