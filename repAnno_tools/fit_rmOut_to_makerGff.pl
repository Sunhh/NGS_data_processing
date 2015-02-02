#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and &stopErr("[Err] perl $0 inFa.rmMsk.out inFa.rmMsk.forMaker.gff\n"); 

my $pl_rm2Gff="/share/app/Annotation/repeatmasker/RepeatMasker/util/rmOutToGFF3.pl"; 
my $gff3_merge="gff3_merge"; 

my $inRmOut = shift; 
my $outGff = shift; 

&exeCmd("perl $pl_rm2Gff $inRmOut > $inRmOut.gff.raw"); 
open F,'<',"$inRmOut.gff.raw" or &stopErr("[Err] Failed to open file [$inRmOut.gff.raw]. $!\n"); 
open O,'>',"tmp.med.gff" or &stopErr("$!\n"); 
my $nn = 0; 
while (<F>) {
	chomp; 
	m/^#/ and next; 
	m/^\>/ and last; 
	my @ta = split(/\t/, $_); 
	$ta[8] =~ m/^Target=/ or &stopErr("[Err] Bad format $ta[8]\n"); 
	$nn++; 
	$ta[8] = "ID=$nn;$ta[8]";
	if ( $ta[8]=~s/(;Target=\S+)\s+0\s+(\d+)/$1 1 $2/ ) {
		&tsmsg("[Rec] Changing Hit_start from 0 to 1 for $ta[8]\n"); 
	}
	print O join("\t", @ta)."\n"; 
}
close O; 
close F; 
&exeCmd("$gff3_merge -l -o $outGff tmp.med.gff"); 
&exeCmd("rm tmp.med.gff"); 


