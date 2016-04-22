#!/usr/bin/perl -w
# 2016-04-21 Strand information is also needed for strand-specific RNA-seq data. 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

!@ARGV and die "perl $0 in.rel_loc in1.agp in2.agp > out.rel_loc\n"; 

my $fn_loc = shift; 
my $fn_agp1 = shift; 
my $fn_agp2 = shift; 

my %c2s_1 = %{ &fileSunhh::load_agpFile( $fn_agp1 ) }; 
my %c2s_2 = %{ &fileSunhh::load_agpFile( $fn_agp2 ) }; 

my $fh = &openFH( $fn_loc, '<' ); 
# head -4 P1toP3_s2.maf.rel_loc
# Cmo_Scf00009    5985    Cma_Scf00001    785     +
# Cmo_Scf00009    5986    Cma_Scf00001    786     +
# Cmo_Scf00009    5987    Cma_Scf00001    787     +
# Cmo_Scf00009    5988    Cma_Scf00001    788     +
my (@b1, @b2, $oStr); 
while (&wantLineC($fh)) {
	my @ta = &splitL("\t", $_); 
	( defined $c2s_1{$ta[0]} and defined $c2s_2{$ta[2]} ) or die "$_\n"; 
	$ta[4] //= '+'; 
	@b1 = $ms_obj->switch_position( 'qry2ref'=> \%c2s_1, 'qryID'=> $ta[0] , 'qryPos' => $ta[1], 'qryStr'=>$ta[4] ); 
	@b2 = $ms_obj->switch_position( 'qry2ref'=> \%c2s_2, 'qryID'=> $ta[2] , 'qryPos' => $ta[3], 'qryStr'=>$ta[4] ); 
	$oStr = ( $b1[0][2] ne $b2[0][2] ) ? '-' : '+' ; 
	print join("\t", @{$b1[0]}[0,1], @{$b2[0]}[0,1], $oStr)."\n"; 
}
close($fh); 

