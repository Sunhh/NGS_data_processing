#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 

!@ARGV and die "perl $0 in_snp.tbl > in_snp.tbl.maf\nPlease note that the geno_col is 3 by default, it is not 2!\n"; 

my $geno_col = 3; 
my %skipColN; 
for (my $i=2; $i<$geno_col; $i++) {
	$skipColN{$i} = 1; 
}

print join("\t", qw/chr pos cnt_all cnt_2Allele MajorBp MajorCnt MinorBp MinorCnt MAF_all MAF_2Allele/)."\n"; 

my $st_obj = SNP_tbl->new('filename'=>shift, 'skipColN'=>\%skipColN); 
&tsmsg("[Msg] Reading SNP\n"); 
$st_obj->readTbl(); 
&tsmsg("[Msg] Counting genotypes\n"); 
$st_obj->cnt_genotype(); 
&tsmsg("[Msg] Output.\n"); 
for (my $i=0; $i<@{ $st_obj->{'chrColV'} }; $i++) {
	my $chrID = $st_obj->{'chrColV'}[$i]; 
	my $pos   = $st_obj->{'posColV'}[$i]; 
	my @aa = @{ $st_obj->{'cnt_alleleTypeBase'}[$i] }; 
	my @vv = @{ $st_obj->{'cnt_alleleCnt'}[$i] }{@aa}; 
	my ($maf1, $maf2); 
	my $cnt_total1 = &mathSunhh::_sum(@vv); 
	my $cnt_total2 = $cnt_total1; 
	if ( @aa == 1 ) {
		push(@aa, 'NA'); 
		push(@vv, '0'); 
		$maf1 = $maf2 = 'NA'; 
		$cnt_total2 = &mathSunhh::_sum(@vv[0,1]); 
	}
	$maf1 //= sprintf("%.4f", $vv[1]/$cnt_total1) * 100; 
	$maf2 //= sprintf("%.4f", $vv[1]/$cnt_total2) * 100; 
	print join("\t", $chrID, $pos, $cnt_total1, $cnt_total2, $aa[0], $vv[0], $aa[1], $vv[1], $maf1, $maf2)."\n"; 
}



