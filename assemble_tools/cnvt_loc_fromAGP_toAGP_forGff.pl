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
	"help!", 
); 

my %glob;
&prepare_glob();  
sub prepare_glob {
	$glob{'help_txt'} = <<HH; 
################################################################################
# perl $0   -old_agp input_old.ctg2scf.agp   -new_agp input_new.ctg2scf.agp   -old_loc in_old_by_oldChrPos.gff3 
#
# -help 
#
# -old_agp        [filename] .AGP format. 
# -new_agp        [filename] .AGP format. 
#                   old_agp and new_agp should share the same contig sets. 
#
# -old_loc        [filename] Gff file. 
#
#
# -new_loc        [filename] This is the output file name. STDOUT if not given. 
#
#
################################################################################
HH

	$glob{'fh_new_loc'} = \*STDOUT; 
	defined $opts{'new_loc'} and $glob{'fh_new_loc'} = &openFH($opts{'new_loc'}, '>'); 
	for my $fn (qw/old_agp new_agp old_loc/) {
		defined $opts{$fn} or do { &tsmsg("[Err]\n"); &tsmsg("[Err] -$fn needed.\n\n"); &LogInforSunhh::usage($glob{'help_txt'}); }; 
		$glob{"fn_$fn"} = $opts{$fn}; 
	}
}# prepare_glob() 

$opts{'help'} and &LogInforSunhh::usage($glob{'help_txt'}); 

%{$glob{'old_c2s'}} = %{ &fileSunhh::load_agpFile( $glob{'fn_old_agp'} ) }; 
%{$glob{'new_c2s'}} = %{ &fileSunhh::load_agpFile( $glob{'fn_new_agp'} ) }; 

%{$glob{'old_s2c'}} = %{ &fileSunhh::reverse_agpHash($glob{'old_c2s'}) }; 

my @aa_loci = &fileSunhh::load_tabFile( $glob{'fn_old_loc'} , 1 ); 

for my $a1 (@aa_loci) {
	# @$a1 == 0 and do { print { $glob{'fh_new_loc'} } "chr\tpos\tstr\n"; next; }; 
	$a1->[0] =~ m!^\s*(#|$)! and do { print { $glob{'fh_new_loc'} } join("\t", @$a1)."\n"; next; }; 
	my ($old_scfID1, $old_scfPos1, $old_scfStr1) = ($a1->[0], $a1->[3], $a1->[6]); 
	my ($old_scfID2, $old_scfPos2, $old_scfStr2) = ($a1->[0], $a1->[4], $a1->[6]); 
	$old_scfStr1 eq '.' and $old_scfStr1 = '+'; 
	$old_scfStr2 eq '.' and $old_scfStr2 = '+'; 
	my @new_scfInf1 = $ms_obj->transfer_position( 'from_ref2qry' => $glob{'old_s2c'}, 'to_qry2ref' => $glob{'new_c2s'}, 'fromLoc' => [$old_scfID1, $old_scfPos1, $old_scfStr1] ); 
	my @new_scfInf2 = $ms_obj->transfer_position( 'from_ref2qry' => $glob{'old_s2c'}, 'to_qry2ref' => $glob{'new_c2s'}, 'fromLoc' => [$old_scfID2, $old_scfPos2, $old_scfStr2] ); 
	$new_scfInf1[2] eq $new_scfInf2[2] or &stopErr("[Err] Bad strand for [@new_scfInf1] [@new_scfInf2]\n"); 
	$new_scfInf1[0] eq $new_scfInf2[0] or &stopErr("[Err] Bad scfID  for [@new_scfInf1] [@new_scfInf2]\n"); 
	$a1->[0] = $new_scfInf1[0]; 
	$a1->[3] = $new_scfInf1[1]; 
	$a1->[4] = $new_scfInf2[1]; 
	$a1->[6] = $new_scfInf1[2]; 
	print { $glob{'fh_new_loc'} } join("\t", @$a1)."\n"; 
}

