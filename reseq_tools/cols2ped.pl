#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 

!@ARGV and die "perl $0 outPref in.cols\n"; 

my $st = SNP_tbl->new(); 

my $oPref = shift; 
my $oPedFile = "$oPref.ped"; 
my $oMapFile = "$oPref.map"; 

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


open OM,'>',"$oMapFile" or die; 
# .map format : chrID, SNP_ID, cM_num(0), position
my (@header, @data, @loci); 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($. == 1) {
		@header = @ta; 
		next; 
	}
	my %cnt_allele; 
	my @tmp_nn; 
	for (my $i=2; $i<@ta; $i++) {
		my @nums = &geno2num($ta[$i]); 
		for my $tn (@nums) {
			$tn == 0 and next; 
			$cnt_allele{$tn} ++; 
		}
		$tmp_nn[$i] = [@nums]; 
	}
	my @srt_allele = sort { $cnt_allele{$b} <=> $cnt_allele{$a} } keys %cnt_allele; 
	if (scalar(@srt_allele) <= 1) {
		# Not only 
		next; 
	} else {
		my %use; 
		for my $tn (@srt_allele[0,1]) {
			$use{$tn} = 1; 
		}
		for ( my $i=2; $i<@tmp_nn; $i++) {
			my $is_0 = 0; 
			for my $tn (@{$tmp_nn[$i]}) {
				defined $use{$tn} or do { $is_0 = 1; last; }; 
			}
			$is_0 == 1 and @{$tmp_nn[$i]} = (0,0); 
			# push(@{$data[$i]}, [@{$tmp_nn[$i]}]); 
			push(@{$data[$i]}, "$tmp_nn[$i][0] $tmp_nn[$i][1]"); 
		}
		push(@loci, [@ta[0,1]]); 
		my $chrID = $ta[0]; $chrID =~ s/^chr(\d+)$/$1/i; $chrID =~ s/^WM97_Chr0*//i; $chrID eq '' and $chrID = 20; 
		print OM join("\t", $chrID, "s${chrID}_$ta[1]", 0, $ta[1])."\n"; 
	}
}
close OM; 
# output .ped 
open OP,'>',"$oPedFile" or die; 
for (my $i=2; $i<@header; $i++) {
	my $pedID = $header[$i]; 
	my $idvID = $header[$i]; 
	my $fatID = 0; 
	my $monID = 0; 
	my $sexID = 1; 
	my $case  = 0; 
	my $genotypeLine = join("\t", @{$data[$i]}); 
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

