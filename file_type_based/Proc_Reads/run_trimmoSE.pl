#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 in_RS.fq out_RS.fq\n"; 

my $inRSFq = shift; 
my $outRSFq = shift; 

my $exe_java = '/share/nas2/xigua/sunhonghe/src/java/jre1.8.0_131/bin/java'; 
my $jar_trim = '/share/nas2/xigua/sunhonghe/src/reads/trimmo/Trimmomatic-0.36/trimmomatic-0.36.jar'; 

my $cmd = ""; 
$cmd .= "$exe_java "; 
$cmd .= " -jar $jar_trim "; 
$cmd .= " SE "; 
$cmd .= " -threads 2 "; 
$cmd .= " $inRSFq $outRSFq "; 
$cmd .= " ILLUMINACLIP:/share/nas2/xigua/sunhonghe/src/reads/trimmo/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 "; 

&exeCmd_1cmd($cmd); 
&exeCmd_1cmd("gzip $outRSFq"); 


