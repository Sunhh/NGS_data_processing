#!/usr/bin/perl 
# This script is designed to count estimates by (per) site . 
# 2015-10-16 Trying to merge some scripts' functions, and clean up the directroy. 
#            Add MAF calculation. 
use strict; 
use warnings; 
use mathSunhh; 
use SNP_tbl; 
my $st_obj = SNP_tbl->new(); 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"startColN:i", # 2 
	"noHeader!", 
); 
$opts{'startColN'} //= 2; 

my $geno_col = $opts{'startColN'}; 

my $help_txt = <<HH; 

Count genotyped (not missing) number, homozygous ratio and heteryzygous ratio per site. 
Also count numbers for MAF. 
Out format : qw/chr pos GenoNum HomoRatio HeteRatio cnt_all cnt_2Allele MajorBp MajorCnt MinorBp MinorCnt MAF_all MAF_2Allele/

perl $0 in.snp_tbl > in.snp_tbl.cntHomHetR

-startColN       [$opts{'startColN'}]

Please note that the geno_col is $opts{'startColN'}

HH
!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

unless ($opts{'noHeader'}) {
	my $l = <>; 
}
print join("\t",qw/chr pos GenoNum HomoRatio HeteRatio cnt_all cnt_2Allele MajorBp MajorCnt MinorBp MinorCnt MAF_all MAF_2Allele/)."\n"; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my $cnt_h = &cnt_in_arr( [@ta[ $opts{'startColN'} .. $#ta ]] ); 
	my $homR = sprintf( "%.4f", $cnt_h->{'homN'}/$cnt_h->{'tot_typed'} * 100 ); 
	my $hetR = sprintf( "%.4f", $cnt_h->{'hetN'}/$cnt_h->{'tot_typed'} * 100 ); 
	print STDOUT join("\t", $ta[0], $ta[1], $cnt_h->{'tot_typed'}, $homR, $hetR, @{$cnt_h}{qw/alleleCnt_all alleleCnt_2 majorBp majorCnt minorBp minorCnt maf_all maf_2/})."\n"; 
}

# chr pos GenoNum HomoRatiio HeteRatio cnt_all cnt_2Allele MajorBp MajorCnt MinorBp MinorCnt MAF_all MAF_2Allele
sub cnt_in_arr {
	# $_[0] : [$geno1, $geno2, $geno3, ... ]
	# $_[1] : Ploidy number for MAF counting. 2 by default. Could be only 1 or 2. 
	# Only count ploidy == 1 or 2 ; 
	$_[1] //= 2; 
	$_[1] == 2 or $_[1] == 1 or &stopErr("[Err] Bad ploidy number [$_[1]] which could be only 1 or 2.\n"); 
	my %cnt; 
	for (@{$_[0]}) {
		$cnt{'tot'} ++; 
		($_ eq 'N' or $_ eq 'n') and do { $cnt{'N'}++; next; }; 
		if ( $_ =~ m/^[ATGC*]$/i or $_ =~ m/\+/ ) {
			$cnt{'homN'}++; $cnt{'tot_typed'}++; 
			$cnt{'allele2cnt'}{$_} += $_[1]; 
		} elsif ( ( my @bb = &SNP_tbl::dna_d2b( &SNP_tbl::dna_b2d($_) ) ) > 0 ) {
			@bb   == 1 and do { $cnt{'homN'}++; $cnt{'tot_typed'}++; $cnt{$bb[0]} += $_[1]; next; }; 
			$_[1] == 1 and do { &tsmsg("[Wrn] Bad genotype [$_] for ploidy=1.\n"); $cnt{'N'}++; next; }; 
			@bb   >  2 and do { &tsmsg("[Wrn] Bad genotype [$_] for ploidy=2.\n"); $cnt{'N'}++; next; }; 
			$cnt{'hetN'}++; $cnt{'tot_typed'}++; 
			for my $tb ( @bb ) { $cnt{'allele2cnt'}{$tb}++; } 
		} else {
			&tsmsg("[Wrn] Weired genotype [$_] is treated as homozygous.\n"); 
			$cnt{'homN'}++; $cnt{'tot_typed'}++; 
			$cnt{'allele2cnt'}{$_} += $_[1]; 
		}
	}# End for (@{$_[0]})

	for (qw/homN hetN N/) {
		$cnt{$_} //= 0; 
	}

	for (qw/tot tot_typed/) {
		$cnt{$_} //= -1; 
	}
	my @aa = sort { $cnt{'allele2cnt'}{$b} <=> $cnt{'allele2cnt'}{$a} || $a cmp $b } keys %{$cnt{'allele2cnt'}}; 
	my @vv = @{ $cnt{'allele2cnt'} }{@aa}; 
	@aa == 0 and do { @aa = ('NA'); @vv = (0) }; 
	if (@aa == 1) {
		push(@aa, 'NA'); 
		push(@vv, 0); 
		$cnt{'maf_all'} = 'NA'; 
		$cnt{'maf_2'}   = 'NA'; 
	}
	$cnt{'majorBp'} = $aa[0]; $cnt{'majorCnt'} = $vv[0]; 
	$cnt{'minorBp'} = $aa[1]; $cnt{'minorCnt'} = $vv[1]; 
	$cnt{'alleleCnt_all'} = $cnt{'tot_typed'} * $_[1]; 
	$cnt{'alleleCnt_all'} == &mathSunhh::_sum( @vv ) or &tsmsg("[Wrn] Bug here! aCnt_all = $cnt{'alleleCnt_all'} for : @{$_[0]}\n"); 
	$cnt{'alleleCnt_all'} == 0 and $cnt{'alleleCnt_all'} = -1; 
	$cnt{'alleleCnt_2'}   = $vv[0] + $vv[1]; $cnt{'alleleCnt_2'} == 0 and $cnt{'alleleCnt_2'} = -1; 
	$cnt{'maf_all'} //= sprintf("%.4f", $vv[1]/$cnt{'alleleCnt_all'}*100); 
	$cnt{'maf_2'}   //= sprintf("%.4f", $vv[1]/$cnt{'alleleCnt_2'}  *100); 

	return (\%cnt); # qw/tot tot_typed N homN hetN alleleCnt_all alleleCnt_2 majorBp majorCnt minorBp minorCnt maf_all maf_2/
}# sub cnt_in_arr() 
# chr pos GenoNum HomoRatiio HeteRatio cnt_all cnt_2Allele MajorBp MajorCnt MinorBp MinorCnt MAF_all MAF_2Allele



