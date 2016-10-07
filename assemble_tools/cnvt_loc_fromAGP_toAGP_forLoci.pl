#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use LogInforSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"old_agp:s", 
	"new_agp:s", 
	"old_loc:s", 
	"new_loc:s", 
	"colN_seqID:s", 
	"colN_seqP:s", 
	"colN_seqStr:s", 
	"help!", 
); 



my %glob;

&prepare_glob();  

%{$glob{'old_c2s'}} = %{ &fileSunhh::load_agpFile( $glob{'fn_old_agp'} ) }; 
%{$glob{'new_c2s'}} = %{ &fileSunhh::load_agpFile( $glob{'fn_new_agp'} ) }; 

%{$glob{'old_s2c'}} = %{ &fileSunhh::reverse_agpHash($glob{'old_c2s'}) }; 

my @aa_loci = &fileSunhh::load_tabFile( $glob{'fn_old_loc'} , 1 ); 

for my $a1 (@aa_loci) {
	# @$a1 == 0 and do { print { $glob{'fh_new_loc'} } "chr\tpos\tstr\n"; next; }; 
	$a1->[0] =~ m!^\s*(#|$)! and do { print { $glob{'fh_new_loc'} } join("\t", @$a1)."\n"; next; }; 
	my ($new_scfID, $new_scfPos, $new_scfStr); 
	my @new_value; 
	for my $cN ( @{$glob{'colN_seqP'}} ) {
		my ($old_scfID, $old_scfPos) = ( @{$a1}[ $glob{'colN_seqID'}[0], $cN ] ); 
		my $old_scfStr = '+'; 
		@{$glob{'colN_seqStr'}} > 0 and $old_scfStr = $a1->[$glob{'colN_seqStr'}[0]]; 
		$old_scfStr eq '.' and $old_scfStr = '+'; 
		my @new_scfInf = $ms_obj->transfer_position( 'from_ref2qry' => $glob{'old_s2c'}, 'to_qry2ref' => $glob{'new_c2s'}, 'fromLoc' => [$old_scfID, $old_scfPos, $old_scfStr] ); 
		$new_scfID  //= $new_scfInf[0]; 
		$new_scfStr //= $new_scfInf[2]; 
		$new_scfID  eq $new_scfInf[0] or &stopErr("[Err] Different new_scfID  for [$old_scfID, $old_scfPos, $old_scfStr]\n"); 
		$new_scfStr eq $new_scfInf[2] or &stopErr("[Err] Different new_strStr for [$old_scfID, $old_scfPos, $old_scfStr]\n"); 
		$new_scfPos = $new_scfInf[1]; 
		$new_value[ $glob{'colN_seqID'}[0] ] //= $new_scfID; 
		@{$glob{'colN_seqStr'}} > 0 and $new_value[ $glob{'colN_seqStr'}[0] ] //= $new_scfStr; 
		$new_value[ $cN ] = $new_scfPos; 
	}
	@{$glob{'colN_seqStr'}} > 0 and $a1->[ $glob{'colN_seqStr'}[0] ] = $new_value[ $glob{'colN_seqStr'}[0] ]; 
	for my $cN ( $glob{'colN_seqID'}[0] , @{$glob{'colN_seqP'}} ) {
		$a1->[$cN] = $new_value[ $cN ]; 
	}
	print { $glob{'fh_new_loc'} } join("\t", @$a1)."\n"; 
}

sub prepare_glob {
	$glob{'colN_seqID'} = '0'; 
	$glob{'colN_seqP'}  = '3,4'; 
	$glob{'colN_seqStr'}= ''; 
	$glob{'help_txt'} = <<HH; 
################################################################################
# perl $0   -old_agp input_old.ctg2scf.agp   -new_agp input_new.ctg2scf.agp   -old_loc in_old_by_oldChrPos.gff3 
#
# -help 
#
# -old_loc        [filename] Input .gff3 file. 
# -new_loc        [filename] Output .gff3 file. print to STDOUT if not given. 
#
# -old_agp        [filename] .AGP format for old_loc
# -new_agp        [filename] .AGP format for new_loc
#                   old_agp and new_agp should share the same contig sets. 
#
# -colN_seqID     ['$glob{'colN_seqID'}']
# -colN_seqP      ['$glob{'colN_seqP'}']
# -colN_seqStr    ['$glob{'colN_seqStr'}']
#
################################################################################
HH

	$opts{'help'} and &LogInforSunhh::usage($glob{'help_txt'}); 

	$glob{'fh_new_loc'} = \*STDOUT; 
	defined $opts{'new_loc'} and $glob{'fh_new_loc'} = &openFH($opts{'new_loc'}, '>'); 
	for my $fn (qw/old_agp new_agp old_loc/) {
		defined $opts{$fn} or do { &tsmsg("[Err]\n"); &tsmsg("[Err] -$fn needed.\n\n"); &LogInforSunhh::usage($glob{'help_txt'}); }; 
		$glob{"fn_$fn"} = $opts{$fn}; 
	}
	for my $tk (qw/colN_seqID colN_seqP colN_seqStr/) {
		defined $opts{$tk} and $glob{$tk} = $opts{$tk}; 
	}
	for my $tk (qw/colN_seqID colN_seqP colN_seqStr/) {
		$glob{$tk} = [ &mathSunhh::_parseCol($glob{$tk}) ]; 
	}
	@{$glob{'colN_seqID'}}  == 1 or &stopErr("[Err] Should set one and only one column number for seqID  [-colN_seqID]\n"); 
	@{$glob{'colN_seqStr'}} <= 1 or &stopErr("[Err] Should set one column number for seqStr at most [-colN_seqStr]\n"); 
	@{$glob{'colN_seqP'}}   >= 1 or &stopErr("[Err] At least set one column number for seqPosition [-colN_seqP]\n"); 
}# prepare_glob() 


