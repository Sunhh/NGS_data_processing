#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 in_RS.fq out_RS.fq\n"; 

my $inRSFq = shift; 
my $outRSFq = shift; 

my $exe_java = '/share/nas2/xigua/sunhonghe/src/java/jre1.8.0_131/bin/java'; 
my $jar_trim = '/share/nas2/xigua/sunhonghe/src/reads/trimmo/Trimmomatic-0.36/trimmomatic-0.36.jar'; 
my $dir_trim = '/data/Sunhh/src/general/trimmomatic/Trimmomatic-0.38/';
my $fn_adp    = "$dir_trim/adapters/TruSeq3-SE.fa";
my $para_trim = " ILLUMINACLIP:${fn_adp}\:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ";

$exe_java = '/usr/lib/jvm/java-11-openjdk-amd64/bin/java';
$jar_trim = '/data/Sunhh/src/general/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar';
$dir_trim = '/data/Sunhh/src/general/trimmomatic/Trimmomatic-0.38/';
$fn_adp    = "$dir_trim/adapters/TruSeq3-SE.fa";
$para_trim = " ILLUMINACLIP:${fn_adp}\:2:30:10 SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 ";

my $cmd = ""; 
$cmd .= "$exe_java "; 
$cmd .= " -jar $jar_trim "; 
$cmd .= " SE "; 
$cmd .= " -threads 2 "; 
$cmd .= " $inRSFq $outRSFq "; 
$cmd .= " $para_trim ";

&exeCmd_1cmd($cmd); 
&exeCmd_1cmd("gzip $outRSFq"); 


