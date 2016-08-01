#!/usr/bin/perl
# 20160729 Add multi-threading 
use strict; 
use warnings; 
use LogInforSunhh; 
use Getopt::Long; 
use Parallel::ForkManager; 
use fileSunhh; 
use fastaSunhh; 
my $fas_obj = fastaSunhh->new(); 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"in_cog:s", # all_orthomcl.ete3.cog
	"in_prot:s", # all_orthomcl.ete3.prot.fa
	"max_geneN:i", # 9999 
	"out_dir:s", # ete3_out/
	"ete3_cmd:s", # 'ete3 build --cpu 25 -w phylomedb4 '
	"cpuN:i", # 1
); 

$opts{'ete3_cmd'}  //= 'ete3 build --cpu 25 -w phylomedb4 '; 
$opts{'max_geneN'} //= 200; 
$opts{'cpuN'}      //= 1; 

my $help_txt = <<HH; 
##################################################################################################
# perl $0 -in_cog all_orthomcl.ete3.cog   -in_prot all_orthomcl.ete3.prot.fa   -out_dir ete3_out/
#
#   This tool is used to generate gene tree within each COG using ETE3 toolkit. 
#
# Reaured : 
# -in_cog        [filename] Format: sp1_seqA1 \\t sp2_seqA2 \\t ... \\n sp1_seqB \\t sp2_seqB \\t ... \\n ...
# -in_prot       [filename] fasta 
# -out_dir       [DirName]  
#
# Optional : 
#
# -max_geneN     [$opts{'max_geneN'}] Maximal number of genes within a OG. 
# -ete3_cmd      ['$opts{'ete3_cmd'}']
#
# -cpuN          [$opts{'cpuN'}] Multi-threading. 
#
# -help 
#
#
####################################################################################################
HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 
defined $opts{'out_dir'} or &LogInforSunhh::usage($help_txt); 
$opts{'out_dir'} =~ s!/+$!!; 

my %glob; 

&tsmsg("[Rec] Loading protein sequence [$opts{'in_prot'}]\n"); 
my %prot_fas = %{ $fas_obj->save_seq_to_hash( 'faFile' => $opts{'in_prot'} ) }; 
for (keys %prot_fas) { chomp($prot_fas{$_}{'seq'}); } 
&tsmsg("[Rec] Loading OG clusters [$opts{'in_cog'}]\n"); 
my @cog_tab = &fileSunhh::load_tabFile( $opts{'in_cog'} ) ; 
-d "$opts{'out_dir'}" or mkdir($opts{'out_dir'}) or &stopErr("[Err] Failed to create out_dir [$opts{'out_dir'}]\n"); 
-d "$opts{'out_dir'}/ete_config/" or mkdir("$opts{'out_dir'}/ete_config/") or &stopErr("[Err] Failed to create out_dir [$opts{'out_dir'}/ete_config/]\n");

my $wrk_dir = &fileSunhh::new_tmp_dir(); 
mkdir($wrk_dir) or &stopErr("[Err] Failed to create dir [$wrk_dir]\n"); 


$glob{'curr_order'} = 0; 
for ( my $i=0; $i<@cog_tab; $i++ ) {
	$glob{'curr_order'} = $i+1; 
	my $curr_geneN = scalar( @{$cog_tab[$i]} ); 
	$curr_geneN > 0 or do { &tsmsg("[Wrn] Skip $glob{'curr_order'} -th COG because no gene found.\n"); next; }; 
	$curr_geneN > $opts{'max_geneN'} and do { &tsmsg("[Wrn] Skip $glob{'curr_order'} -th COG because of too many genes [$curr_geneN]\n"); next; }; 

	### Prepare COG 
	$glob{'curr_tab'} = [ $cog_tab[$i] ]; 
	&prepare_ete3(); 
}#End for ( my $i=0; $i<@cog_tab; $i++ ) 

### Run ete 
$glob{'MAX_PROCESSES'} = $opts{'cpuN'}; 
$glob{'pm'} = new Parallel::ForkManager($glob{'MAX_PROCESSES'}); 
for my $sep ( @{$glob{'sepPref'}} ) {

	my $has_td = ( defined $glob{'tree_dir'} ) ? 1 : 0; 
	$has_td == 1 and $glob{'pid'} = $glob{'pm'}->start and next; 

	my %tGlob = %$sep; 
	my $pref = $tGlob{'curr_pref'}; 
	&exeCmd_1cmd("$opts{'ete3_cmd'} -a ${pref}.prot.fa -o ${pref}_ete3O/", 0) and &stopErr("[Err] Failed to run ete3 for $tGlob{'curr_order'} -th COG\n"); 
	opendir DD, "${pref}_ete3O/" or &stopErr("[Err] Failed to opendir [${pref}_ete3O/]\n"); 
	for my $ff (readdir(DD)) {
		$ff =~ m/(^\.|^db$|^tasks$)/i and next; 
		if ( !(defined $glob{'tree_dir'}) and -d "${pref}_ete3O/$ff" ) {
			$glob{'tree_dir'} = $ff; 
			-d "$opts{'out_dir'}/$glob{'tree_dir'}/" or mkdir("$opts{'out_dir'}/$glob{'tree_dir'}/") or &stopErr("[Err] Failed to create dir [$opts{'out_dir'}/$glob{'tree_dir'}/]"); 
		}
		if ( $ff =~ m!^ete_build\.cfg$!i ) {
			&fileSunhh::_move( "${pref}_ete3O/$ff", "$opts{'out_dir'}/ete_config/$tGlob{'curr_order'}.$ff" ); 
		} elsif ( $ff eq $glob{'tree_dir'} ) {
			opendir(D2, "${pref}_ete3O/$ff/") or &stopErr("[Err] Faield to opendir [${pref}_ete3O/$ff]\n"); 
			for my $f2 ( readdir(D2) ) {
				$f2 =~ m!^\.! and next; 
				if ( $f2 =~ m!^$tGlob{'curr_order'}\.! ) {
					&fileSunhh::_move( "${pref}_ete3O/$ff/$f2", "$opts{'out_dir'}/$glob{'tree_dir'}/$f2" ); 
				} elsif ( $f2 =~ m!^(runid|command_lines|commands.log)!i ) {
					&fileSunhh::_move( "${pref}_ete3O/$ff/$f2", "$opts{'out_dir'}/$glob{'tree_dir'}/$tGlob{'curr_order'}.$f2" );
				} else {
					&stopErr("[Err] Undefined file [${pref}_ete3O/$ff/$f2] in ete output dir [${pref}_ete3O/$ff/]\n"); 
				}
			}
			closedir(D2); 
			for my $fn (@{$tGlob{'curr_toRM'}}) {
				unlink($fn); 
			}
		} else {
			&stopErr("[Err] Undefined [${pref}_ete3O/$ff] in ete output dir [${pref}_ete3O/]\n"); 
		}
	}# End for my $ff (readdir(DD))
	closedir(DD); 
	&fileSunhh::_rmtree("${pref}_ete3O/"); 

	$has_td == 1 and $glob{'pm'}->finish; 

}#End for my $pref () 
$glob{'pm'}->wait_all_children; 

&fileSunhh::_rmtree($wrk_dir); 

&tsmsg("[Rec] All done [$0]\n"); 

sub prepare_ete3 {
	# A proper $glob{'curr_order'} is required. 
	### To be used : curr_order , curr_tab 
	for my $tk (qw/curr_pref curr_fh_oProt curr_toRM/) {
		delete($glob{$tk}); 
	}
	$glob{'curr_pref'} = "$wrk_dir/$glob{'curr_order'}"; 
	&tsmsg("[Msg]   Prepared data for [$glob{'curr_pref'}]\n"); 
	$glob{'curr_fh_oProt'} = &openFH("$glob{'curr_pref'}.prot.fa", '>'); 
	my %h; 
	for my $ar ( @{$glob{'curr_tab'}} ) {
		for my $id (@$ar) {
			defined $h{$id} and next; 
			defined $prot_fas{$id} or &stopErr("[Err] No prot-seq found for ID [$id]\n"); 
			print {$glob{'curr_fh_oProt'}} ">$prot_fas{$id}{'key'}\n$prot_fas{$id}{'seq'}\n"; 
		}
	}
	delete($glob{'curr_tab'}); 
	close( $glob{'curr_fh_oProt'} ); delete( $glob{'curr_fh_oProt'} ); 
	push(@{$glob{'curr_toRM'}}, "$glob{'curr_pref'}.prot.fa"); 
	push(@{$glob{'file_toRM'}}, @{$glob{'curr_toRM'}}); 
	my %t_glob = %glob; 
	delete($t_glob{'file_toRM'}); 
	delete($t_glob{'sepPref'}); 
	push( @{$glob{'sepPref'}}, \%t_glob ); 
}#sub prepare_ete3()


