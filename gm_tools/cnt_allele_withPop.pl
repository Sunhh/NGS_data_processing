#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use SNP_tbl; 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"vcf_tab:s", # Input file coming from "vcf_to_tab"; The first three columns are #CHROM , POS and REF . 
	"P1_colN:i", # 3
	"P2_colN:i", # 4
); 

$opts{'P1_colN'} //= 3; 
$opts{'P2_colN'} //= 4; 

my $help_txt = <<HH; 

perl $0 -vcf_tab Itay_Galap_132offs_2parents_geno.tbl 

-help

-P1_colN       [$opts{'P1_colN'}]
-P2_colN       [$opts{'P2_colN'}]

HH

defined $opts{'vcf_tab'} or &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %glob; 

$glob{'fh_inTab'} = &openFH($opts{'vcf_tab'}, '<'); 
$glob{'need_header'} = 1; 

while (readline($glob{'fh_inTab'})) {
	chomp($_); 
	my @ta = &splitL("\t", $_); 
	if ($ta[0] =~ m!^(#CHROM|chr|chromosome|chrID|CHROM)$!i) {
		# This is header. 
		$glob{'tab_header'} //= [ @ta ]; 
		next; 
	}
	unless ( defined $glob{'offs_colN'} ) {
		for (my $i=3; $i<@ta; $i++) {
			$i == $opts{'P1_colN'} and next; 
			$i == $opts{'P2_colN'} and next; 
			push(@{$glob{'offs_colN'}}, $i); 
		}
	}
	# Get P1 allele : 
	my @p1_al = &SNP_tbl::tab_allele( $ta[$opts{'P1_colN'}] ); 
	my @p2_al = &SNP_tbl::tab_allele( $ta[$opts{'P2_colN'}] ); 
	my %pp_al_class = %{ &SNP_tbl::tab_class_PP_al( \@p1_al, \@p2_al ) }; 
	
	# Class offsprings
	my %cnt; 
	my %genoCnt; 
	for my $i (@{$glob{'offs_colN'}}) {
		$cnt{'total_N'} ++; 
		my @of_al = &SNP_tbl::tab_allele( $ta[$i] ); 
		if ($#of_al == 0) {
			$genoCnt{"$of_al[0][0]/$of_al[0][0]"} ++; 
		} else {
			$genoCnt{"$of_al[0][0]/$of_al[1][0]"} ++; 
		}
		my $of_class = &SNP_tbl::tab_class_off_al( \%pp_al_class, \@of_al ); 
		if ( $of_class eq 'miss' ) {
			$cnt{'miss_N'} ++; 
		} elsif ( $of_class eq 'homo_P1_parent' ) {
			$cnt{'homo_P1'} ++; 
		} elsif ( $of_class eq 'homo_P2_parent' ) {
			$cnt{'homo_P2'} ++; 
		} elsif ( $of_class eq 'hete_both_parent' ) {
			$cnt{'hete_PP'} ++; 
		} elsif ( $of_class =~ m!_non_parent$! ) {
			$cnt{'non_PP'} ++; 
		} elsif ( $of_class =~ m!_bad_parent$! ) {
			$cnt{'bad_PP'} ++; 
		} else {
			$cnt{'other'}{ $of_class } ++; 
		}
	}
	$cnt{'other'}{'genoCnt'} = join(':', map { "|$_|=$genoCnt{$_}" } sort keys %genoCnt); 
	if ( $glob{'need_header'} ) {
		$glob{'tab_header'} //= [ 'chr', 'pos', 'ref' ]; 
		$glob{'tab_header'}[$opts{'P1_colN'}] //= 'P1'; $glob{'tab_header'}[$opts{'P1_colN'}] eq '' and $glob{'tab_header'}[$opts{'P1_colN'}] = 'P1'; 
		$glob{'tab_header'}[$opts{'P2_colN'}] //= 'P2'; $glob{'tab_header'}[$opts{'P2_colN'}] eq '' and $glob{'tab_header'}[$opts{'P2_colN'}] = 'P2'; 
		print STDOUT join("\t", @{$glob{'tab_header'}}[0,1,2, $opts{'P1_colN'}, $opts{'P2_colN'}], qw/total miss homo_P1 homo_P2 hete_PP non_PP bad_PP others/)."\n"; 
		$glob{'need_header'} = 0; 
	}
	for my $tk (qw/total_N miss_N homo_P1 homo_P2 hete_PP non_PP bad_PP/) {
		$cnt{$tk} //= 0; 
	}
	$cnt{'other'} //= { 'NA'=>'NA' }; 
	print STDOUT join("\t", 
	  @ta[0,1,2,$opts{'P1_colN'},$opts{'P2_colN'}], 
	  @cnt{qw/total_N miss_N homo_P1 homo_P2 hete_PP non_PP bad_PP/}, 
	  join(";;", map { "$_=$cnt{'other'}{$_}" } sort keys %{$cnt{'other'}})
	)."\n"; 
}
close($glob{'fh_inTab'}); 

