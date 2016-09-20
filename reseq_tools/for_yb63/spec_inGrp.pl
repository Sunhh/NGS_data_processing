#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_tbl:s", 
	"list_gIn:s", "list_gOut:s", 
	"missR_gIn:f", "missR_gOut:f", 
	"major_ratio_gIn:f", 
	"max_ratio_gOut:f", 
); 

my $help_txt = <<HH; 
perl $0 -in_tbl in_snp.tbl -list_gIn colN_list_of_inner_grp   -list_gOut colN_list_of_outer_grp

-missR_gIn         [0] missing ratio allowed in inner group 
-missR_gOut        [0] missing ratio allowed in outer group 
-major_ratio_gIn   [1] major genotype allowed in inner group
-max_ratio_gOut    [0] genotype ratio allowed in outer group

HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'in_tbl'} or &LogInforSunhh::usage($help_txt);
defined $opts{'list_gIn'} or &LogInforSunhh::usage($help_txt);
defined $opts{'list_gOut'} or &LogInforSunhh::usage($help_txt);


$opts{'missR_gIn'}  //= 0; 
$opts{'missR_gOut'} //= 0; 
$opts{'major_ratio_gIn'} //= 1; 
$opts{'max_ratio_gOut'}  //= 0; 

$opts{'major_ratio_gIn'} > 0.5 or die "-major_ratio_gIn is supposed to bigger than 0.5\n"; 

my @colN_gIn  = map { $_->[0] } &fileSunhh::load_tabFile( $opts{'list_gIn'} ); 
my @colN_gOut = map { $_->[0] } &fileSunhh::load_tabFile( $opts{'list_gOut'} ); 
my %glob; 
$glob{'chk_gIn' } = scalar(@colN_gIn ); 
$glob{'chk_gOut'} = scalar(@colN_gOut); 

my $fh = &openFH($opts{'in_tbl'}, '<'); 
my @hh; 
{
	my $h = <$fh>; chomp($h); @hh = split(/\t/, $h); 
	print STDOUT "$h\n"; 
}
SNP: 
while (<$fh>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my %cnt; 
	$cnt{'cntN_in'} = 0; 
	$cnt{'cntN_out'} = 0; 
	for my $i (@colN_gIn) {
		$ta[$i] = uc($ta[$i]); 
		$ta[$i] = join('', sort split(//, $ta[$i])); 
		$ta[$i] eq 'N' and do { $cnt{'cntN_in'}++; next; }; 
		$ta[$i] =~ m/^[ATGC]+$/ or die "|$ta[$i]|\n"; 
		$cnt{'g2n_in'}{$ta[$i]} ++; 
	}
	$glob{'chk_gIn'} - $cnt{'cntN_in'} > 0 or next SNP; 
	$cnt{'cntN_in'} <= $opts{'missR_gIn'} * $glob{'chk_gIn'} or next SNP; 
	for my $geno1 (sort { $cnt{'g2n_in'}{$b} <=> $cnt{'g2n_in'}{$a} || $a cmp $b } keys %{$cnt{'g2n_in'}}) {
		push(@{$cnt{'g2n_in_arr'}}, [ $geno1, $cnt{'g2n_in'}{$geno1} ]); 
	}
	my $spec_geno = $cnt{'g2n_in_arr'}[0][0]; 
	$cnt{'g2n_in_arr'}[0][1] >= $opts{'major_ratio_gIn'} * ($glob{'chk_gIn'} - $cnt{'cntN_in'}) or next SNP; 
	for my $i (@colN_gOut) {
		$ta[$i] = uc($ta[$i]); 
		$ta[$i] = join('', sort split(//, $ta[$i])); 
		$ta[$i] eq 'N' and do { $cnt{'cntN_out'}++; next; }; 
		$ta[$i] =~ m/^[ATGC]+$/ or die "|$ta[$i]|\n"; 
		$cnt{'g2n_out'}{$ta[$i]} ++; 
	}
	$glob{'chk_gOut'} - $cnt{'cntN_out'} > 0 or next SNP; 
	$cnt{'cntN_out'} <= $opts{'missR_gOut'} * $glob{'chk_gOut'} or next SNP; 
	if ( defined $cnt{'g2n_out'}{ $spec_geno } ) {
		$cnt{'g2n_out'}{ $spec_geno } <= $opts{'max_ratio_gOut'} * ( $glob{'chk_gOut'} - $cnt{'cntN_out'} ) or next SNP; 
	}
	print STDOUT "$_\n"; 
}
close ($fh); 


