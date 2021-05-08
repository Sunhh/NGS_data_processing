#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, "help!", 
	"P1ID:s", 
	"P2ID:s", 
	"B1ID:s", 
	"minGQ_P12:i", 
	"minDP_P12:i", 
	"minDP_P1:i", 
	"minDP_P2:i", 
	"minDP_B1:i", 
); 

$opts{'minGQ_P12'} //= -1; 
$opts{'minDP_P12'} //= -1; 
$opts{'minDP_P1'}  //= -1; 
$opts{'minDP_P2'}  //= -1; 

$opts{'minDP_B1'}  //= -1; 

my $help_txt = <<HH; 
################################################################################
# perl $0 -P1ID SYNLvWangTuo   -P2ID BWJMoLvTuo   -B1ID F2BsaR lmyPM2019bsa_filtV_PASS.vcf > lmyPM2019bsa_filtV_PASS.vcf.bsa_SNPidx
#
# -minGQ_P12      [$opts{'minGQ_P12'}]
# -minDP_P12      [$opts{'minDP_P12'}]
# -minDP_B1       [$opts{'minDP_B1'}]
#
################################################################################
HH

-t and !@ARGV and &LogInforSunhh::usage($help_txt); 


my %inf; 
LINE: 
while (<>) {
	chomp; 
	if (m!^\s*#!) {
		# m!^\s*##! and do { print STDOUT "$_\n"; next LINE; }; 
		m!^\s*##! and do { next LINE; }; 
		if ( m!^#CHROM\t! ) {
			my @ta=split(/\t/, $_); 
			for (my $i=9; $i<@ta; $i++) {
				for my $id (qw/P1ID P2ID B1ID/) {
					defined $opts{$id} and $ta[$i] eq $opts{$id} and $inf{'id_cN'}{$id} = $i; 
				}
			}
			for my $id (qw/P1ID P2ID B1ID/) {
				defined $opts{$id} or next LINE; 
				defined $inf{'id_cN'}{$id} or &stopErr("[Err] Failed to find ID for [$id] [$opts{$id}]\n"); 
			}
			print STDOUT join("\t", qw/CHROM POS REF ALT/, "P1_$opts{'P1ID'}", "P2_$opts{'P2ID'}", $opts{'B1ID'}, "$opts{'B1ID'}:SNPi")."\n"; 
			next LINE; 
		} else {
			&stopErr("[Err] Bad line: $_\n"); 
		}
	}
	my @ta=split(/\t/, $_); 
	my $is_good = 1; 
	my $vR_P1 = &val_genotype( $ta[8], $ta[ $inf{'id_cN'}{'P1ID'} ] ); $vR_P1->{'e'} == 0 or do { &tsmsg("[Wrn] Skip bad line [e=$vR_P1->{'e'}]: $_\n"); next LINE; }; 
	my $vR_P2 = &val_genotype( $ta[8], $ta[ $inf{'id_cN'}{'P2ID'} ] ); $vR_P2->{'e'} == 0 or do { &tsmsg("[Wrn] Skip bad line [e=$vR_P2->{'e'}]: $_\n"); next LINE; }; 

	$vR_P1->{'GT'} =~ m!^(\d+)/\1$!         or $is_good = 0; $is_good == 0 and next LINE;  
	$vR_P1->{'GQ'} >= $opts{'minGQ_P12'}    or $is_good = 0; $is_good == 0 and next LINE;  
	$vR_P1->{'DP'} >= $opts{'minDP_P12'}    or $is_good = 0; $is_good == 0 and next LINE;  
	$vR_P2->{'GT'} =~ m!^(\d+)/\1$!         or $is_good = 0; $is_good == 0 and next LINE;  
	$vR_P2->{'GQ'} >= $opts{'minGQ_P12'}    or $is_good = 0; $is_good == 0 and next LINE;  
	$vR_P2->{'DP'} >= $opts{'minDP_P12'}    or $is_good = 0; $is_good == 0 and next LINE;  
	$vR_P1->{'GT'} eq $vR_P2->{'GT'} and next LINE; 

	my $vR_B1 = &val_genotype( $ta[8], $ta[ $inf{'id_cN'}{'B1ID'} ] ); $vR_B1->{'e'} == 0 or do { &tsmsg("[Wrn] Skip bad line [e=$vR_B1->{'e'}]: $_\n"); next LINE; };
	$vR_B1->{'e'} == 2 and do { &tsmsg("[Wrn] Skip bad line: $_\n"); next LINE; }; 
	$vR_B1->{'GT'} eq './.' and next LINE; 
	$vR_B1->{'DP'} >= $opts{'minDP_B1'} or next LINE; 
	my ($allele_1) = ( $vR_P1->{'GT'} =~ m!^(\d+)\/\1$! ); 
	my ($allele_2) = ( $vR_P2->{'GT'} =~ m!^(\d+)\/\1$! ); 
	my @B1_cntAll  = split(/,/, $vR_B1->{'AD'}); 
	my $B1_cnt1 = $B1_cntAll[$allele_1]; 
	my $B1_cnt2 = $B1_cntAll[$allele_2]; 
	$B1_cnt1 + $B1_cnt2 == 0 and next; 
	my $B1_snpIdx = $B1_cnt1/($B1_cnt1+$B1_cnt2); 
	print STDOUT join("\t", @ta[0,1,3,4], $vR_P1->{'GT'}, $vR_P2->{'GT'}, "$vR_B1->{'GT'}:$vR_B1->{'AD'}", $B1_snpIdx)."\n"; 
}

sub val_genotype {
	my ($fmt, $geno) = @_; 
	my %back; 
	$back{'e'} = 0; 
	my @a_f = split(/:/, $fmt); 
	my @a_g = split(/:/, $geno); 
	# scalar(@a_f) == scalar(@a_g) or &stopErr("[err] bad fmt_wi_geno [$fmt] [$geno]\n"); 
	scalar(@a_f) == scalar(@a_g) or do { $back{'e'}=2; return (\%back);  }; 
	for (my $i=0; $i<@a_f; $i++) {
		$back{$a_f[$i]} = $a_g[$i]; 
	}
	return(\%back); 
}


