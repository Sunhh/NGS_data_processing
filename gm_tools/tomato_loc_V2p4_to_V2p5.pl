#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use LogInforSunhh; 

!@ARGV and die "perl $0 in_old.agp in_new.agp loci_list\n"; 

my $fn_oldAGP = shift; 
my $fn_newAGP = shift; 
my $fn_loci = shift; 

my %old_c2s = %{ &fileSunhh::load_agpFile( $fn_oldAGP ) }; 
my %new_c2s = %{ &fileSunhh::load_agpFile( $fn_newAGP ) }; 

my %old_s2c = %{ &fileSunhh::reverse_agpHash(\%old_c2s) }; 
for my $sID (keys %old_s2c) {
	@{$old_s2c{ $sID }} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$old_s2c{ $sID }}; 
}

my @aa_loci = &fileSunhh::load_tabFile( $fn_loci , 1 ); 
for my $a1 (@aa_loci) {
	@$a1 == 0 and do { print "chr\tpos\tstr\n"; next; }; 
	$a1->[0] =~ m!^\s*#! and do { print join("\t", @$a1, qw/chr pos str/)."\n"; next; }; 
	my ($old_scfID, $old_scfPos) = ($a1->[0], $a1->[1]); 
	my @new_scfInf = $ms_obj->transfer_position( 'from_ref2qry' => \%old_s2c, 'to_qry2ref' => \%new_c2s, 'fromLoc' => [$old_scfID, $old_scfPos] ); 
	print join("\t", $old_scfID, $old_scfPos, $new_scfInf[0], $new_scfInf[1], $new_scfInf[2])."\n"; 
}

