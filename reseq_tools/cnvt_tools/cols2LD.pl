#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

my $pl_dir = '/data/Sunhh/watermelon/source_reseq/new_source/04.LD/hap_LD'; 
my $pl_rmSame = "$pl_dir/" . 'rm_same_site_hete2N.pl'; 
my $pl_mkHap = "$pl_dir/" . 'tbl2hap.pl'; 
my $pl_binLD = "$pl_dir/" . 'get_the_LD_decay_file.pl'; 

!@ARGV and die "perl $0 in.snp\n"; 

my $snpF = shift; 

&exeCmd_1cmd("perl $pl_rmSame $snpF > $snpF.var"); 
&exeCmd_1cmd("perl $pl_mkHap $snpF.var $snpF.ped $snpF.info"); 
&exeCmd_1cmd("java -Xmx20960M -jar /data/Sunhh/src/Evolution/haploview/Haploview4.1.jar -nogui -minMAF 0.05 -hwcutoff 0.001 -dprime -log $snpF.log -out $snpF -pedfile $snpF.ped -info $snpF.info"); 
&exeCmd_1cmd("perl $pl_binLD $snpF.LD 1000 $snpF.LD_bin1k"); 

&tsmsg("[Rec] Done.\n"); 

