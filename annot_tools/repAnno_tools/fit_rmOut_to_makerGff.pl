#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"pl_rm2Gff:s", 
	"exe_gff3_merge:s", 
	"help!", 
); 

$opts{'pl_rm2Gff'} //= '/Data/Sunhh/src/annotation/repeatmasker/RepeatMasker/util/rmOutToGFF3.pl'; # NERCV server; 
$opts{'exe_gff3_merge'} //= 'gff3_merge'; 

my $htxt = <<HH; 

perl $0 inFa.rmMsk.out inFa.rmMsk.forMaker.gff > aaaa

 -pl_rm2Gff       [Path to rmOutToGFF3.pl] $opts{'pl_rm2Gff'}
 -exe_gff3_merge  [gff3_merge] From maker. 

HH

!@ARGV and &LogInforSunhh::usage($htxt); 
$opts{'help'} and &LogInforSunhh::usage($htxt); 


my $inRmOut = shift; 
my $outGff = shift; 



my $pl_rm2Gff = $opts{'pl_rm2Gff'}; 
my $gff3_merge = $opts{'exe_gff3_merge'}; 

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


