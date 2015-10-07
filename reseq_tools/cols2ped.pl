#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	'help!', 
	'max2allele!', 
	'startColN:i', 
	'noHeader!', 
); 
$opts{'startColN'} //= 2; 

sub usage {
	print STDERR <<HH; 

perl $0 outPref in.cols\n
-help 
-max2allele     Only accept sites with two alleles if given; 
-startColN      [2] 0-indexed number of the first column for sample base. 
-noHeader       Given if input .cols file has no header. 

HH
	exit 1; 
}

!@ARGV and &usage(); 
$opts{'help'} and &usage(); 

my $st = SNP_tbl->new(); 

my $oPref = shift; 
my $oPedFile = "$oPref.ped"; 
my $oMapFile = "$oPref.map"; 
my $oInfFile = "$oPref.info"; 

my %dblist; 
{
my @aa = (
[qw/W A T/], 
[qw/S C G/],
[qw/M A C/],
[qw/K G T/],
[qw/R A G/],
[qw/Y C T/]
); 
for my $tr (@aa) {
	my @bb = @$tr; 
	$dblist{$bb[0]} = [$bb[1], $bb[2]]; 
}
}


# .ped file format : 
#   pedID indvID FatherID MotherID sexTag caseTag genotypes(A=1, C=2, G=3, T=4, Miss=0)
my %allele2num = qw(
A 1
C 2
G 3
T 4
N 0
); 
my %num2allele = qw(
1 A
2 C
3 G
4 T
0 N
); 


open OM,'>',"$oMapFile" or die; 
open OI,'>',"$oInfFile" or die; 
# .map format : chrID, SNP_ID, cM_num(0), position
my (@header, @data); 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1 and !($opts{'noHeader'})) {
		@header = @ta; 
		next; 
	}
	my %cnt_allele; 
	my @tmp_nn; 
	for (my $i=$opts{'startColN'}; $i<@ta; $i++) {
		my @nums = &geno2num($ta[$i]); 
		for my $tn (@nums) {
			$tn == 0 and next; 
			$cnt_allele{$tn} ++; 
		}
		$tmp_nn[$i] = [@nums]; 
	}
	if ($opts{'max2allele'}) {
		scalar(keys %cnt_allele) == 2 or next; 
	} else {
		scalar(keys %cnt_allele) <= 1 and next; 
	}
	my @srt_allele = sort { $cnt_allele{$b} <=> $cnt_allele{$a} } keys %cnt_allele; 
	my %use; 
	for my $tn (@srt_allele[0,1]) {
		$use{$tn} = 1; 
	}
	for ( my $i=$opts{'startColN'}; $i<@tmp_nn; $i++) {
		my $is_0 = 0; 
		for my $tn (@{$tmp_nn[$i]}) {
			defined $use{$tn} or do { $is_0 = 1; last; }; 
		}
		$is_0 == 1 and @{$tmp_nn[$i]} = (0,0); 
		push(@{$data[$i]}, "$tmp_nn[$i][0] $tmp_nn[$i][1]"); 
	}
	my $chrID = $ta[0]; $chrID =~ s/^chr(\d+)$/$1/i; $chrID =~ s/^WM97_Chr0*//i; $chrID eq '' and $chrID = 200; 
	print OM join("\t", $chrID, "s${chrID}_$ta[1]", 0, $ta[1])."\n"; 
	my $af_0 = sprintf("%.4f", $cnt_allele{$srt_allele[0]}/( $cnt_allele{$srt_allele[0]} + $cnt_allele{$srt_allele[1]}) ); 
	my $af_1 = sprintf("%.4f", $cnt_allele{$srt_allele[1]}/( $cnt_allele{$srt_allele[0]} + $cnt_allele{$srt_allele[1]}) ); 
	print OI join("\t", $#{$data[$opts{'startColN'}]}+1, $ta[1], $num2allele{ $srt_allele[0] }, $af_0, $num2allele{ $srt_allele[1] }, $af_1)."\n"; 
}
close OM; 
close OI; 
# output .ped 
open OP,'>',"$oPedFile" or die; 
for (my $i=0; $i+$opts{'startColN'}<@data; $i++) {
	my $pedID = ($opts{'noHeader'}) ? $i+1 : $header[$i+$opts{'startColN'}]; 
	my $idvID = ($opts{'noHeader'}) ? $i+1 : $header[$i+$opts{'startColN'}]; 
	my $fatID = 0; 
	my $monID = 0; 
	my $sexID = 0; # 0/1? 
	my $case  = 0; 
	my $genotypeLine = join("\t", @{$data[$i+$opts{'startColN'}]}); 
	print OP join("\t", $pedID, $idvID, $fatID, $monID, $sexID, $case, $genotypeLine)."\n"; 
}
close OP; 

sub geno2num {
	my $tb = shift; 
	$tb = $st->SingleChar($tb); 
	if (defined $allele2num{$tb}) {
		return ($allele2num{$tb}, $allele2num{$tb}); 
	} elsif ( defined $dblist{$tb} ) {
		my @b1 = @{$dblist{ $tb }}; 
		return ($allele2num{$b1[0]}, $allele2num{$b1[1]}); 
	} else {
		return (0, 0); 
	}
}

