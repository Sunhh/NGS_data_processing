#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use SNP_tbl; 
-t and !@ARGV and die "perl $0 in_tbl.snp\n"; 
my $fn = shift; 

my $st = SNP_tbl->new('filename'=>$fn); 

&tsmsg("[Rec] Reading $fn\n"); 
$st->readTbl(); 
&tsmsg("[Rec] First SingleCharData()\n"); 
$st->SingleCharData('onlyATGC'=>0, 'maxAlleleN'=>0); # Remove indels in the table. 
&tsmsg("[Rec] First count genotypes.\n"); 
$st->cnt_genotype(); 
&tsmsg("[Rec] Recording Values.\n"); 
print STDOUT join("\t", qw/Chrom Pos cnt_lmiss cnt_alleleTypeN alleleTypeCnt cnt_homo cnt_hete rmHete_cnt_lmiss rmHete_cnt_alleleTypeN rmHete_alleleTypeCnt rmHete_cnt_homo rmHete_MAF/)."\n"; 
# 0       Chrom
# 1       Pos
# 2       cnt_lmiss
# 3       cnt_alleleTypeN
# 4       alleleTypeCnt
# 5       cnt_homo
# 6       cnt_hete
# 7       rmHete_cnt_lmiss
# 8       rmHete_cnt_alleleTypeN
# 9       rmHete_alleleTypeCnt
# 10      rmHete_cnt_homo
# 11      rmHete_MAF
my @out_arr; 
for (my $i=0; $i<@{$st->{'data_arr'}}; $i++) {
	push(@{$out_arr[$i]}, 
	  $st->{'chrColV'}[$i], # For Chrom
	  $st->{'posColV'}[$i], # Pos
	  $st->{'cnt_lmiss'}[$i], # cnt_lmiss, number of individuals missing in current site. 
	  $st->{'cnt_alleleTypeN'}[$i], # cnt_alleleTypeN, number of allele types in current site. 
	  join(';;', map { "${_}_$st->{'cnt_alleleCnt'}[$i]{$_}"; } @{ $st->{'cnt_alleleTypeBase'}[$i] } ), 
	  $st->{'cnt_homo'}[$i], # cnt_homo, number of homozygous individuals in current site. 
	  $st->{'cnt_hete'}[$i], # cnt_hete, number of heterozygous individuals in current site. 
	); 
}
&tsmsg("[Rec] Convert heterozygous genotypes to N.\n"); 
$st->SingleCharData('onlyATGC'=>1); # Remove heterozygous in the table. 
&tsmsg("[Rec] Second count genotype.\n"); 
$st->cnt_genotype(); 
&tsmsg("[Rec] Count MAF\n"); 
$st->cnt_maf(); 
&tsmsg("[Rec] Record values again.\n"); 
for (my $i=0; $i<@{$st->{'data_arr'}}; $i++) {
	push(@{$out_arr[$i]}, 
	  $st->{'cnt_lmiss'}[$i], # rmHete_cnt_lmiss 
	  $st->{'cnt_alleleTypeN'}[$i], # rmHete_cnt_alleleTypeN 
	  join(';;', map { "${_}_$st->{'cnt_alleleCnt'}[$i]{$_}" } @{ $st->{'cnt_alleleTypeBase'}[$i] } ), 
	  $st->{'cnt_homo'}[$i], 
	  $st->{'cnt_maf'}[$i]
	); 
}
&tsmsg("[Rec] Output values' table\n"); 
for my $oa (@out_arr) {
	print STDOUT join("\t", @$oa)."\n"; 
}
&tsmsg("[Rec] Finished.\n"); 
