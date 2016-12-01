#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use SNP_tbl; 

@ARGV >= 3 or die "perl $0 in.snp out_hap.ped out_hap.info\n"; 

my $inTbl = shift; 
my $outTbl = shift; 
my $outInfo = shift; 

my $st = SNP_tbl->new('filename' => $inTbl ); 
$st->readTbl(); 
$st->tbl2hap( 'ofile'=>$outTbl, 'oinfofile'=>$outInfo ); 

for my $a ( 1 .. 11) {
	my $chrID = sprintf("WM97_Chr%02d", $a); 
	my $file = "$chrID/total.$chrID.snp"; 
	
}

