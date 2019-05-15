#!/usr/bin/perl
use strict; 
use warnings; 


# [Sunhh@bioinfor01 201905]$ head -5 WM_SNP_byWM97V1chr_segPos.fa.topbChrV1.bn6
# WMSNP141_Left100        WM97pbV1_Chr01  100.000 101     0       0       1       101     81975   82075   9.14e-47        187     101     36935898        plus
# WMSNP142_1_Left100      WM97pbV1_Chr01  100.000 101     0       0       1       101     2198389 2198489 9.14e-47        187     101     36935898        plus
# WMSNP143_Left100        WM97pbV1_Chr01  100.000 101     0       0       1       101     6102235 6102335 9.14e-47        187     101     36935898        plus
# WMSNP144_Left100        WM97pbV1_Chr01  100.000 101     0       0       1       101     8083296 8083396 9.14e-47        187     101     36935898        plus
# WMSNP145_Left100        WM97pbV1_Chr01  100.000 101     0       0       1       101     26344185        26344085        9.14e-47        187     101     36935898        minus

my %ppp; 
while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	$ta[2] >= 100 or next; 
	my $str = '+'; 
	if ($ta[8] > $ta[9]) {
		# @ta[8,9] = @ta[9,8]; 
		$str = '-'; 
	}
	my ($snpID); 
	if ($ta[0] =~ m!^(\S+)_Left100$!) {
		$snpID = $1; 
		$ta[7] == $ta[12] or next; 
		$ta[7]-$ta[6]+1 >= 0.8*$ta[12] or next; 
		if (defined $ppp{$snpID}{'left'}) {
			$ppp{$snpID}{'left'}[1] eq 'R' and next; 
			$ppp{$snpID}{'left'}[1] = 'R'; 
		} else {
			$ppp{$snpID}{'left'} = [$ta[1], $ta[9], $str]; 
		}
	} elsif ($ta[0] =~ m!^(\S+)_Right100$!) {
		$snpID = $1; 
		$ta[6] == 1 or next; 
		$ta[7]-$ta[6]+1 >= 0.8*$ta[12] or next; 
		if (defined $ppp{$snpID}{'right'}) {
			$ppp{$snpID}{'right'}[1] = 'R'; 
		} else {
			$ppp{$snpID}{'right'} = [$ta[1], $ta[8], $str]; 
		}
	} else {
		die "$_\n"; 
	}
}
print STDOUT join("\t", qw/SNP_ID Final_chr Final_pos Final_strand left_chr left_pos left_strand right_chr right_pos right_strand/)."\n"; 
for my $snpID (sort keys %ppp) {
	$ppp{$snpID}{'left'}  //= ['NA','NA','NA']; 
	$ppp{$snpID}{'right'} //= ['NA','NA','NA']; 
	my $isGood = 1; 
	for my $lr (qw/left right/) {
		$ppp{$snpID}{$lr}[1] eq 'R' and $isGood = 0; 
		$ppp{$snpID}{$lr}[0] eq 'NA' and $isGood = 0; 
	}
	if ($isGood == 1) {
		for (my $i=0; $i<3; $i++) {
			$ppp{$snpID}{'left'}[$i] eq $ppp{$snpID}{'right'}[$i] or $isGood = 0; 
		}
	}
	my ($oChr, $oPos, $oStr) = qw/NA NA NA/; 
	if ($isGood == 1) {
		$oChr = $ppp{$snpID}{'left'}[0]; 
		$oPos = $ppp{$snpID}{'left'}[1]; 
		$oStr = $ppp{$snpID}{'left'}[2]; 
	}
	print STDOUT join("\t", $snpID, $oChr, $oPos, $oStr, @{$ppp{$snpID}{'left'}}[0..2], @{$ppp{$snpID}{'right'}}[0..2])."\n"; 
}



