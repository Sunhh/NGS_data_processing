#!/usr/bin/perl
use strict; 
use warnings; 
use fileSunhh; 
use SNP_tbl; 

!@ARGV and die "perl $0 out_pref parent_LA483_UC204B_132offs.comm.snp\n"; 
my $pref = shift; 


# ==> parent_LA483_UC204B_132offs.comm.snp <==
# chr     pos     base    LA483   UC-204B Itay19  Itay2   Itay55  Itay75  RIL-1-G RIL-1-R RIL-11-G        RIL-11R RIL-12-G        RIL-12-R        RIL-13-G        RIL-13-
# SL2.40ch00      549274  T       G/G     T/T     T/T     T/T     T/T     G/G     T/T     T/T     T/T     T/G     ./.     T/T     ./.     T/T     T/T     T/T     T/T
# SL2.40ch00      551416  C       T/T     C/C     C/C     C/C     C/C     ./.     C/C     C/C     ./.     ./.     C/C     ./.     ./.     C/C     ./.     C/C     C/C
# SL2.40ch00      1146669 C       T/T     C/C     C/C     C/C     C/C     ./.     C/C     ./.     C/C     C/C     ./.     ./.     C/C     C/C     C/C     C/C     C/C


my $hd = <>; 
{
	chomp($hd); 
	my @ta=split(/\t/, $hd); 
	my $b1 = "marker\tposition(bp)"; 
	my $b2 = join("\t", @ta[5 .. $#ta]); 
	$hd = join("\t", $b1, $b2)."\n"; 
}

my %have; 

while (<>) {
	chomp; 
	my @ta=split(/\t/, $_); 
	my $b1 = "$ta[0]_$ta[1]\t$ta[1]"; 
	my $fn = "${pref}$ta[0].csv"; 
	unless (defined $have{$fn}) {
		$have{$fn} = 1; 
		&fileSunhh::write2file($fn, $hd, '>' ); 
	}
	my @p1_al = &SNP_tbl::tab_allele( $ta[3] ); 
	my @p2_al = &SNP_tbl::tab_allele( $ta[4] ); 
	( $#p1_al == 0 and $#p2_al == 0 ) or next; 
	my %pp_class = %{ &SNP_tbl::tab_class_PP_al( \@p1_al, \@p2_al ) }; 
	my $is_bad = 0; 
	for my $tc ( @ta[5..$#ta] ) {
		my @of_al = &SNP_tbl::tab_allele( $tc ); 
		my $of_class = &SNP_tbl::tab_class_off_al( \%pp_class, \@of_al ); 
		if ( $of_class eq 'miss' ) {
			$tc = '-'; 
		} elsif ( $of_class eq 'homo_P1_parent' ) {
			$tc = 'a'; 
		} elsif ( $of_class eq 'homo_P2_parent' ) {
			$tc = 'b'; 
		} elsif ( $of_class eq 'hete_both_parent' ) {
			$tc = 'h'; 
		} else {
			$is_bad = 1; 
			last; 
		}
	}
	$is_bad == 1 and next; 
	&fileSunhh::write2file($fn, join("\t", $b1, @ta[5 .. $#ta])."\n", '>>'); 
}

