#!/usr/bin/perl
use strict; 
use fileSunhh; 
use LogInforSunhh; 

# [Sunhh@Penguin GeneticMap_V2_RILseq]$ head -5 WM97_vFix1.scf2chr.list LIST_20110305 | more
# ==> WM97_vFix1.scf2chr.list <==
# Chr_id  Scaffold_id     Oriented        Start   End
# chr1    scaffold1481    F       1       1224693
# chr1    scaffold9256    N       1225694 1315378
# chr1    scaffold27      F       1316379 4288155
# chr1    scaffold1911    R       4289156 5260582
#
# ==> LIST_20110305 <==
# Chr1    scaffold1481    Y       F       1       1224693
# Chr1    scaffold27      Y       F       1225694 4197470
# Chr1    scaffold619     Y       R       4198471 4893635
# Chr1    scaffold1911    Y       N       4894636 5866062
# Chr1    scaffold1546    Y       N       5867063 6904819

use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"scf2chr_lis1:s", 
	"scf2chr_lis2:s", 
); 


my %str2num; 
%str2num = qw(
 1       1
 -1      -1
 f       1
 r       -1
 plus    1
 minus   -1
 n       0
); 

my $lfh1 = &openFH($opts{'scf2chr_lis1'}, '<'); 
my $lfh2 = &openFH($opts{'scf2chr_lis2'}, '<'); 
my %s2c_1; # {scf} = [chrID, strand, start, end, scf_length]; 
my %c2s_1; # {chr} = {loci}{scfID} = [ strand, start, end, scf_length]; {scfID}{list} = [ scfID1, scfID2, ...]
while (<$lfh1>) {
	chomp; m/^\s*$|^\s*#/ and next; 
	my @ta = split(/\t/, $_); 
	my ($chrID, $scfID, $strand, $cS, $cE) = @ta; 
	$strand = lc($strand); 
	defined $str2num{$strand} or die "strand=$strand unknown\n"; 
	$strand = $str2num{$strand}; 
	$cS > $cE and ($cS, $cE) = ($cE, $cS); 
	push(@{$c2s_1{$chrID}{'loci'}{$scfID}}, [ $strand, $cS, $cE, $cE-$cS+1 ]); 
	push(@{$c2s_1{$chrID}{'list'}}, $scfID); 
	$s2c_1{$scfID} = [ $chrID, $strand, $cS, $cE, $cE-$cS+1 ]; 
}
my %s2c_2; 
my %c2s_2; 
while (<$lfh2>) {
	chomp; m/^\s*$|^\s*#/ and next; 
	my @ta = split(/\t/, $_); 
	my ($chrID, $scfID, $strand, $cS, $cE) = @ta; 
	$strand = lc($strand); 
	defined $str2num{$strand} or die "strand=$strand unknown\n"; 
	$strand = $str2num{$strand}; 
	$cS > $cE and ($cS, $cE) = ($cE, $cS); 
	push(@{$c2s_2{$chrID}{'loci'}{$scfID}}, [ $strand, $cS, $cE, $cE-$cS+1 ]); 
	push(@{$c2s_2{$chrID}{'list'}}, $scfID); 
	$s2c_2{$scfID} = [ $chrID, $strand, $cS, $cE, $cE-$cS+1 ]; 
}
close($lfh1); 
close($lfh2); 

# Sort scfID by loci
for my $chrID (keys %c2s_1) {
	@{$c2s_1{$chrID}{'list'}} = sort { $c2s_1{$chrID}{'loci'}{$a}[1] <=>$c2s_1{$chrID}{'loci'}{$b}[1]  } @{$c2s_1{$chrID}{'list'}}; 
}
for my $chrID (keys %c2s_2) {
	@{$c2s_2{$chrID}{'list'}} = sort { $c2s_2{$chrID}{'loci'}{$a}[1] <=>$c2s_2{$chrID}{'loci'}{$b}[1]  } @{$c2s_2{$chrID}{'list'}}; 
}

while (<>) {
	chomp; 
	if ( m/^\s*($|#)/ ) {
		print STDOUT "$_\n"; 
		next; 
	}
	my @ta = split(/\t/, $_); 
	if ($ta[0] =~ m/^(chr|ChromID)$/i) {
		print STDOUT "$_\n"; 
		next; 
	}
	my ($chrID, $chrPos) = @ta[0,1]; 
	defined $c2s_1{$chrID} or &stopErr("[Err] Undefined chrID=[$chrID] in list 1\n"); 
	my ($scfID, $scfPos, $scfStr) = &chrPosInScf($chrID, $chrPos, \%c2s_1, 1); 
	my ($v2_chrID, $v2_chrPos, $v2_chrStr) = &scfPosInChr($scfID, $scfPos, $scfStr); 
	if ( $v2_chrStr == -1 ) {
		for my $tb (@ta[2 .. $#ta]) {
			$tb =~ tr/acgturykmbvdhACGTURYKMBVDHwWsSnN/tgcaayrmkvbhdTGCAAYRMKVBHDwWsSnN/; 
			$tb = reverse($tb); 
			# $tb =~ s/^([ATGCN]+)(\*|\+)([ATGCN]+)$/$3$2$1/; 
			$tb = "rev:$tb"; 
		}
	} elsif ( $v2_chrStr == 1 ) {
	} else {
		&stopErr("[Err] Bad v2_chrStr [$v2_chrStr]\n"); 
	}
	print STDOUT join("\t", @ta)."\n"; 
}

sub chrPosInScf {
	my ($chrID, $chrPos, $c2s_href, $chrStr) = @_; 
	defined $chrStr or $chrStr = 1; 
	defined $c2s_href->{$chrID} or &stopErr("[Err] chrID [$chrID] not defined in c2s_href\n"); 
	for my $scfID (@{$c2s_href->{$chrID}{'list'}}) {
		$c2s_href->{$chrID}{'loci'}{$scfID}[1] > $chrPos and next; 
		$c2s_href->{$chrID}{'loci'}{$scfID}[3] < $chrPos and next; 
		my $back_scfID = $scfID; 
		my ( $strand, $cS, $cE, $scfLen ) = @{ $c2s_href->{$chrID}{'loci'}{$scfID} };
		my $back_scfPos = ''; 
		my $back_scfStr = $chrStr; 
		$back_scfStr =~ m/^0+$/ and $back_scfStr = 1; 
		if ( $strand == -1 ) {
			$back_scfPos = $cE - $chrPos + 1; 
			$back_scfStr *= -1; 
		} elsif ( $strand =~ m/^0|1$/ ) {
			$back_scfPos = $chrPos - $cS + 1; 
		} else {
			&stopErr("[Err] unknown strand number [$strand]\n"); 
		}
		
		return ( $back_scfID, $back_scfPos, $back_scfStr ); 
	}
	return ['', '', '']; 
}

sub scfPosInChr {
	my ( $scfID, $scfPos, $s2c_href, $scfStr ) = @_; 
	defined $scfStr or $scfStr = 1; 
	defined $s2c_href->{$scfID} or &stopErr("[Err] [$scfID] not defined in s2c_href\n"); 
	my ($chrID, $strand, $cS, $cE, $scfLen) = @{ $s2c_href->{$scfID} }; 
	my ($back_chrID, $back_chrPos, $back_chrStr) = ('', '', ''); 
	$back_chrID = $chrID; 
	$back_chrStr = $scfStr; 
	$back_chrStr =~ m/^0+$/ and $back_chrStr = 1; 
	if ( $strand == -1 ) {
		$back_chrPos = $cE - $scfPos + 1; 
		$back_chrStr *= -1; 
	} elsif ( $strand =~ m/^0|1$/ ) {
		$back_chrPos = $cS + $scfPos - 1; 
	} else {
		&stopErr("[Err] unknown scfStr number [$scfStr]\n"); 
	}
	return ( $back_chrID, $back_chrPos, $back_chrStr ); 
}


