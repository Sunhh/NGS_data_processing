#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

!@ARGV and die "\nperl $0 P1toP3_s2.maf.rel_loc.chr_newP3 P1Genom_V1p2.prot_chr.gff3.bed_fit_YZ_wiCDS\n\n"; 

my $file_1to1 = shift; 
my $file_bed  = shift; 

my %loc_1to1; 
&load_loc_1to1( $file_1to1, \%loc_1to1 ); 

my $fh_bed = &openFH( $file_bed, '<' ); 
while (&fileSunhh::wantLineC( $fh_bed )) {
	my ( $oChrID, $oChrS, $oChrE, $oGeneID, $oCDSLen, $oStr, $oCDSLoc ) = split(/\t/, $_); 
	my ( $nChrID1, $nChrS, $nStr1) = &get_newLoc( \%loc_1to1, $oChrID, $oChrS, $oStr ); 
	my ( $nChrID2, $nChrE, $nStr2) = &get_newLoc( \%loc_1to1, $oChrID, $oChrE, $oStr ); 
	my @oCDSLoc_arr = map { [ split(/,/, $_) ] } split(/;/, $oCDSLoc); 
	my @nCDSLoc_arr; 
	if ( $nChrID1 eq $nChrID2 and $nStr1 eq $nStr2 and $nChrID1 ne '' ) {
		for my $t1 (@oCDSLoc_arr) {
			my ( $t2_ID_1, $t2_S, $t2_str_1 ) = &get_newLoc( \%loc_1to1, $oChrID, $t1->[0], $oStr ); 
			my ( $t2_ID_2, $t2_E, $t2_str_2 ) = &get_newLoc( \%loc_1to1, $oChrID, $t1->[1], $oStr ); 
			$t2_S eq '' and $t2_S = 'NA'; 
			$t2_E eq '' and $t2_E = 'NA'; 
			$t2_ID_1 ne '' and $t2_ID_1 ne $nChrID1 and $t2_S = 'NA'; 
			$t2_ID_2 ne '' and $t2_ID_2 ne $nChrID1 and $t2_E = 'NA'; 
			$t2_str_1 ne '' and $t2_str_1 ne $nStr1 and $t2_S = 'NA'; 
			$t2_str_2 ne '' and $t2_str_2 ne $nStr1 and $t2_S = 'NA'; 
			$nStr1 ne $oStr and ($t2_S, $t2_E) = ( $t2_E , $t2_S ); 
			push( @nCDSLoc_arr, [$t2_S, $t2_E] ); 
		}
		$nStr1 ne $oStr and ($nChrS, $nChrE) = ($nChrE, $nChrS); 

	} else {
			$nChrID1=$nChrID2='NA'; 
			$nChrS = $nChrE = 'NA'; 
			$nStr1 = $nStr2 = 'NA'; 
			for my $t1 (@oCDSLoc_arr) {
				push( @nCDSLoc_arr, ['NA', 'NA'] ); 
			}
	}

	print STDOUT join("\t", $_, $nChrID1, $nChrS, $nChrE, $nStr1, join(';', map { "$_->[0],$_->[1]" } @nCDSLoc_arr))."\n"; 
}
close($fh_bed); 

sub load_loc_1to1 {
	my ($fn, $hr) = @_;
	$hr //= {};
	my $fh = &openFH( $fn, '<' );
	&tsmsg("[Msg] Loading loc_1to1 [$fn]\n");
	my %c;
	$c{'log_section'} = { 'cntN_base' => 0, 'cntN_step' => 1e6 };
	while (&fileSunhh::wantLineC( $fh )) {
		&fileSunhh::log_section( $. , $c{'log_section'} ) and &tsmsg("[Msg] Reading $. line.\n");
		my ( $oID, $oP, $nID, $nP, $nStr ) = split(/\t/, $_);
		push( @{$hr->{$oID}{$oP}}, [ $nID, $nP, $nStr ] );
	}
	close($fh);
	&tsmsg("[Msg] Loaded [$fn]\n");
	return(\$hr);
}# load_loc_1to1()

sub get_newLoc {
	my ($hr, $oID, $oP, $oStr) = @_; # input Old coordinates.
	$oStr =~ m/^[+-]$/ or &stopErr("[Err] Bad strand input of get_newLoc( $hr, $oID, $oP, $oStr ).\n");
	my ($nID, $nP, $nStr) = ('', '', '');
	my $flank_len = 50;
	defined $hr->{$oID} or return($nID, $nP);
	if (defined $hr->{$oID}{$oP}) {
		($nID, $nP, $nStr) = @{ $hr->{$oID}{$oP}[0] };
	} else {
		for (my $i=1; $i<=$flank_len; $i++) {
			my $tP_1 = $oP-$i;
			defined $hr->{$oID}{$tP_1} and do { ($nID, $nP, $nStr) = @{ $hr->{$oID}{$tP_1}[0] }; last; };
			my $tP_2 = $oP+$i;
			defined $hr->{$oID}{$tP_2} and do { ($nID, $nP, $nStr) = @{ $hr->{$oID}{$tP_2}[0] }; last; };
		}
	}
	$oStr eq '-' and $nStr =~ tr/+-/-+/;
	return($nID, $nP, $nStr);
}# get_newLoc()



