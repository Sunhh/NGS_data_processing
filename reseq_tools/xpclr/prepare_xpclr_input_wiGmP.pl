#!/usr/bin/perl
use strict; 
use fileSunhh; 
use LogInforSunhh; 
use mathSunhh; 
use wm97Sunhh; 
use SNP_tbl; 

!@ARGV and die "perl $0   out_pref   apple.snp_addGmP   idv_list_1_objPop   idv_list_2_refPop\n"; 

my $opref  = shift; 
my $fn_snp = shift; 
my $fn_lis1 = shift; 
my $fn_lis2 = shift;
my $fn_chrID2num = shift; 
$fn_chrID2num //= ''; 


my %glob; 
my $gmCn = 2; 


if ( defined $fn_chrID2num and $fn_chrID2num ne '' ) {
	my $fh = &openFH( $fn_chrID2num , '<' ); 
	while (&wantLineC($fh)) {
		my @ta = &splitL( "\t", $_ ); 
		( defined $ta[0] and defined $ta[1] ) or next; 
		$glob{'chr_id2num'}{$ta[0]} = $ta[1]; 
	}
	close($fh); 
}

my %fho; 
my %fno; 

my %indv_1 = &indv_list($fn_lis1); 
my %indv_2 = &indv_list($fn_lis2); 

my @h; 
my $fh_snp = &openFH( $fn_snp, '<' ); 
{ my $a=&wantLineC($fh_snp); @h=split(/\t/, $a); } 

{
	my %nn; 
	my @new_h; 
	for (my $i=0; $i<@h; $i++) {

		defined $indv_1{'has'}{$h[$i]} and push( @{$indv_1{'goodIdx'}}, $i ); 
		defined $indv_2{'has'}{$h[$i]} and push( @{$indv_2{'goodIdx'}}, $i ); 
		defined $indv_1{'has'}{$h[$i]} and defined $indv_2{'has'}{$h[$i]} and &stopErr("[Err] Indv [$h[$i]] exists in both groups.\n"); 
	}
	@{$indv_1{'goodSample'}} = @h[ @{$indv_1{'goodIdx'}} ]; 
	@{$indv_2{'goodSample'}} = @h[ @{$indv_2{'goodIdx'}} ]; 
}

my %cc = ( "cntN_base"=>0, "cntN_step"=>1e4 ); 

my %locs; 
my %sites; 
my %prev; 
while (&wantLineC($fh_snp)) { 
	&fileSunhh::log_section($. , \%cc) and &tsmsg("[Msg] $. line.\n"); 
	my @ta=split(/\t/, $_); 
#	defined $phy2gm_P{$ta[0]} or next; 
#	defined $phy2gm_P{$ta[0]}{$ta[1]} or next; 

	my @ta_1 = @ta[ @{$indv_1{'goodIdx'}} ]; 
	my @ta_2 = @ta[ @{$indv_2{'goodIdx'}} ]; 

	my ( %al, %al_1, %al_2 ); 
	my %geno; 
	my $is_bad = 0; 
	for (my $i=0; $i<@ta_1; $i++ ) { 
		if ($ta_1[$i] =~ m/^([ATGCN])$/) { 
			$al_1{$1} += 2; 
			$al{$1}   += 2; 
		} elsif ($ta_1[$i] =~ m/^([ATGC])([ATGC])$/) { 
			$al_1{$1} ++; $al_1{$2} ++; 
			$al{$1}   ++; $al{$2}   ++; 
		} elsif ( $ta_1[$i] =~ m/^[ATGC]{3,}$/ ) { 
			$is_bad = 1; 
			last; 
		} else { 
			my @bb = &SNP_tbl::dna_d2b( $ta_1[$i] ); 
			if ( @bb == 2 ) {
				$al_1{$bb[0]} ++; $al_1{$bb[1]} ++; 
				$al{$bb[0]}   ++; $al{$bb[1]}   ++; 
			} else {
				if ( !(defined $glob{'bad_geno'}{$ta_1[$i]}) ) {
					$glob{'bad_geno'}{$ta_1[$i]} = 1; 
					&tsmsg("[Wrn] Skip site with bad genotype [$ta_1[$i]]\n"); 
				}
				$is_bad = 1; 
				last; 
			}
		}
	}
	delete $al_1{'N'}; 
	scalar(keys %al_1) > 0 or $is_bad = 1; 
	$is_bad == 1 and next; 
	for (my $i=0; $i<@ta_2; $i++ ) { 
		if ($ta_2[$i] =~ m/^([ATGCN])$/) { 
			$al_2{$1} += 2; 
			$al{$1}   += 2; 
		} elsif ($ta_2[$i] =~ m/^([ATGC])([ATGC])$/) { 
			$al_2{$2} ++; $al_2{$2} ++; 
			$al{$2}   ++; $al{$2}   ++; 
		} elsif ( $ta_2[$i] =~ m/^[ATGC]{3,}$/ ) { 
			$is_bad = 1; 
			last; 
		} else {
			my @bb &SNP_tbl::dna_d2b( $ta_2[$i] ); 
			if ( @bb == 2 ) {
				$al_2{$bb[0]} ++; $al_2{$bb[1]} ++; 
				$al{$bb[0]}   ++; $al{$bb[1]}   ++; 
			} else {
				if ( !(defined $glob{'bad_geno'}{$ta_2[$i]}) ) {
					$glob{'bad_geno'}{$ta_2[$i]} = 1; 
					&tsmsg("[Wrn] Skip site with bad genotype [$ta_2[$i]]\n"); 
				}
				$is_bad = 1; 
				last; 
			}
		}
	}
	delete $al_2{'N'}; 
	scalar(keys %al_2) > 0 or $is_bad = 1; 
	$is_bad == 1 and next; 
	delete $al{'N'}; 
	scalar( keys %al ) == 2 or next; 
	my @aa = sort { $al{$b} <=> $al{$a} || $a cmp $b } keys %al; 
	$geno{$aa[0]}         = '1 1'; 
	$geno{"$aa[0]$aa[0]"} = '1 1'; 
	$geno{$aa[1]}         = '0 0'; 
	$geno{"$aa[1]$aa[1]"} = '0 0'; 
	$geno{"$aa[0]$aa[1]"} = '1 0'; 
	$geno{"$aa[1]$aa[0]"} = '1 0'; 
	$geno{"N"}            = '9 9'; 

	
	my $chrID = $ta[0]; 
	my $chrN ; 
	if ( defined $glob{'chr_id2num'}{$chrID} ) {
		$chrN = $glob{'chr_id2num'}{$chrID}; 
	} else {
		$chrN = &wm97Sunhh::chrID_to_number( $chrID , 'WM97_Chr'); 
		$chrN =~ m!^\d+$! or &stopErr("[Err] Failed to convert chrID [$chrID] to number [$chrN]\n"); 
		defined $glob{'chr_num2id'}{$chrN} and &stopErr("[Err] Repeat chrN [$chrN] for differnt chrID [$glob{'chr_num2id'}{$chrN} $chrID]\n"); 
		$glob{'chr_id2num'}{$chrID} = $chrN; 
		$glob{'chr_num2id'}{$chrN}  = $chrID; 
	}

	$chrN =~ m!^\d+$! or &stopErr("[Err] Bad chrID [$chrN] from [$ta[0]]\n"); 

	my $gmP = $ta[$gmCn]; 
	if (defined $prev{'gmID'} and $prev{'gmID'} eq $ta[0]) {
		$prev{'gmP'} < $gmP or next; 
	}
	$prev{'gmID'} = $ta[0]; 
	$prev{'gmP'}  = $gmP; 

	#unless ( defined $fho{$ta[0]} ) {
	#	$fho{$ta[0]}{'geno1'} = &openFH( "${opref}.$ta[0]_g1.geno", '>' ); 
	#	$fho{$ta[0]}{'geno2'} = &openFH( "${opref}.$ta[0]_g2.geno", '>' ); 
	#	$fho{$ta[0]}{'SNP'}   = &openFH( "${opref}.$ta[0].snp",     '>' ); 
	#}
	#
	#
	#print {$fho{$ta[0]}{'geno1'}} join(' ', map { $geno{$_} } @ta_1)."\n"; 
	#print {$fho{$ta[0]}{'geno2'}} join(' ', map { $geno{$_} } @ta_2)."\n"; 
	#
	#my $mrkID = "$ta[0]_$ta[1]"; 
#	print {$fho{$ta[0]}{'SNP'}} join( "\t", $mrkID, $chrN, $phy2gm_P{$ta[0]}{$ta[1]}, $ta[1], $aa[0], $aa[1] )."\n"; 
	#print {$fho{$ta[0]}{'SNP'}} join( "\t", $mrkID, $chrN, $gmP                     , $ta[1], $aa[0], $aa[1] )."\n"; 

	unless ( defined $fno{$ta[0]} ) {
		$fno{$ta[0]}{'geno1'} = "${opref}.$ta[0]_g1.geno"; 
		$fno{$ta[0]}{'geno2'} = "${opref}.$ta[0]_g2.geno"; 
		$fno{$ta[0]}{'SNP'}   = "${opref}.$ta[0].snp"; 
		&fileSunhh::write2file( $fno{$ta[0]}{'geno1'}, '', '>' ); 
		&fileSunhh::write2file( $fno{$ta[0]}{'geno2'}, '', '>' ); 
		&fileSunhh::write2file( $fno{$ta[0]}{'SNP'}  , '', '>' ); 
	}

	&fileSunhh::write2file( $fno{$ta[0]}{'geno1'}, join(' ', map { $geno{$_} } @ta_1)."\n", '>>' ); 
	&fileSunhh::write2file( $fno{$ta[0]}{'geno2'}, join(' ', map { $geno{$_} } @ta_2)."\n", '>>' ); 

	my $mrkID = "$ta[0]_$ta[1]"; 
	&fileSunhh::write2file( $fno{$ta[0]}{'SNP'}, join( "\t", $mrkID, $chrN, $gmP                     , $ta[1], $aa[0], $aa[1] )."\n", '>>' ); 

}
#for my $chrID ( keys %fho ) {
#	for my $k2 ( qw/geno1 geno2 SNP/ ) {
#		close( $fho{$chrID}{$k2} ); 
#	}
#}

&tsmsg("[Rec] Done. $0\n"); 

#my $fh_oSite = &openFH("${opref}.sites", '>'); 
#my $fh_oLoci = &openFH("${opref}.locs",  '>'); 
sub indv_list {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		my @ta=&splitL("\t", $_); 
		defined $back{'has'}{$ta[0]} and next; 
		push(@{$back{'arr'}}, $ta[0]); 
		$back{'has'}{$ta[0]} = $#{$back{'arr'}}; 
	}
	close ($fh); 
	return(%back); 
}


sub load_gmP {
	my $fn = shift; 
	my %back; 
	my $fh = &openFH($fn, '<'); 
	while (&wantLineC($fh)) {
		# chr \\t pos \\t cM 
		my @ta = &splitL("\t", $_); 
		$ta[0] eq 'chr' and next; 
		$back{$ta[0]}{$ta[1]} = $ta[2]; 
	}
	close($fh); 
	return(%back); 
}


