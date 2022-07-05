#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 in_R1.fq in_R2.fq out_Prefix\n"; 

my $inR1Fq = shift; 
my $inR2Fq = shift; 
my $outPref = shift; 
my $outR1PFq = "${outPref}_pTr_R1.fq"; 
my $outR1SFq = "${outPref}_sTr_R1.fq"; 
my $outR2PFq = "${outPref}_pTr_R2.fq"; 
my $outR2SFq = "${outPref}_sTr_R2.fq"; 

my $exe_java = '/share/nas2/xigua/sunhonghe/src/java/jre1.8.0_131/bin/java'; 
my $jar_trim = '/share/nas2/xigua/sunhonghe/src/reads/trimmo/Trimmomatic-0.36/trimmomatic-0.36.jar'; 
my $dir_trim = '/Data/Sunhh/src/reads/trimmomatic/Trimmomatic-0.36/'; 
my $fn_adp   = "$dir_trim/adapters/TruSeq3-PE-2.fa"; 
my $para_trim = " ILLUMINACLIP:${fn_adp}\:2:30:10:1:TRUE SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 "; 

$exe_java = '/usr/java/jre1.8.0_144/bin/java'; 
$jar_trim = '/Data/Sunhh/src/reads/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar'; 
$dir_trim = '/Data/Sunhh/src/reads/trimmomatic/Trimmomatic-0.36/'; 
$fn_adp   = "$dir_trim/adapters/TruSeq3-PE-2.fa"; 
$para_trim = " ILLUMINACLIP:${fn_adp}\:2:30:10:1:TRUE SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:40 "; 


my $cmd = ""; 
$cmd .= "$exe_java "; 
$cmd .= " -jar $jar_trim "; 
$cmd .= " PE "; 
$cmd .= " -threads 2 "; 
$cmd .= " $inR1Fq $inR2Fq $outR1PFq $outR1SFq $outR2PFq $outR2SFq "; 
$cmd .= " $para_trim "; 

&exeCmd_1cmd($cmd); 
&exeCmd_1cmd("gzip $outR1PFq $outR1SFq $outR2PFq $outR2SFq"); 


