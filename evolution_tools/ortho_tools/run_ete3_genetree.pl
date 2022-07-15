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
	"in_cds:s", # all_orthomcl.ete3.cds.fa
	"thres_nt_switch_para:s", # --nt-switch-threshold 0.9 
	"max_geneN:i", # 9999 
	"out_dir:s", # ete3_out/
	"ete3_cmd:s", # 'ete3 build --cpu 2 -w standard_fasttree --clearall '
	"redo_list:s", # 
	"cpuN:i", # 1
); 

# $opts{'ete3_cmd'}  //= 'ete3 build --cpu 2 -w phylomedb4 --clearall '; 
$opts{'ete3_cmd'}  //= 'ete3 build --cpu 2 -w standard_fasttree --clearall '; 
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
# -redo_list     [filename] A list of group IDs which need to be done. 
#
# -cpuN          [$opts{'cpuN'}] Multi-threading. 
#
# -in_cds        [filename] fasta
# -thres_nt_switch_para ['--nt-switch-threshold 0.9']
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
if (defined $opts{'redo_list'}) {
	my $fh = &openFH($opts{'redo_list'}, '<'); 
	while (&wantLineC($fh)) {
		my @ta=&splitL("\t", $_); 
		$glob{'redo_ID'}{$ta[0]} //= $.; 
	}
	close($fh); 
}

&tsmsg("[Rec] Loading protein sequence [$opts{'in_prot'}]\n"); 
my %prot_fas = %{ $fas_obj->save_seq_to_hash( 'faFile' => $opts{'in_prot'} ) }; 
for (keys %prot_fas) { chomp($prot_fas{$_}{'seq'}); } 
my %cds_fas; 
if ( defined $opts{'in_cds'} ) {
	&tsmsg("[Rec] Loading CDS sequence [$opts{'in_cds'}]\n"); 
	%cds_fas = %{ $fas_obj->save_seq_to_hash( 'faFile' => $opts{'in_cds'} ) }; 
	for (keys %cds_fas) { chomp($cds_fas{$_}{'seq'}); }
	$opts{'thres_nt_switch_para'} //= "--nt-switch-threshold 0.9"; 
}
&tsmsg("[Rec] Loading OG clusters [$opts{'in_cog'}]\n"); 
my @cog_tab = &fileSunhh::load_tabFile( $opts{'in_cog'} ) ; 
-d "$opts{'out_dir'}" or mkdir($opts{'out_dir'}) or &stopErr("[Err] Failed to create out_dir [$opts{'out_dir'}]\n"); 
-d "$opts{'out_dir'}/ete_config/" or mkdir("$opts{'out_dir'}/ete_config/") or &stopErr("[Err] Failed to create out_dir [$opts{'out_dir'}/ete_config/]\n");

my $wrk_dir = &fileSunhh::new_tmp_dir(); 
mkdir($wrk_dir) or &stopErr("[Err] Failed to create dir [$wrk_dir]\n"); 


for (my $i=0; $i<@cog_tab; $i++) {
	push(@{$glob{'cog_idx2geneN'}}, [$i, $#{$cog_tab[$i]}+1]); 
}
for my $tr ( sort { $a->[1] <=> $b->[1] || $a->[0] <=> $b->[0] } @{$glob{'cog_idx2geneN'}} ) {
	my $i = $tr->[0]; 
	$glob{'curr_order'} = $i + 1; 

	if (defined $glob{'redo_ID'} and (keys %{$glob{'redo_ID'}})>0) {
		defined $glob{'redo_ID'}{ $glob{'curr_order'} } or next; 
	}

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
	&tsmsg(join("","[Msg] Files [", join(" ", @{$tGlob{'curr_toRM'}}),"] will be removed after pref [$pref] processed\n")); 
	if ( defined $opts{'in_cds'} ) {
		&exeCmd_1cmd("$opts{'ete3_cmd'} -a ${pref}.prot.fa -o ${pref}_ete3O/ -n ${pref}.cds.fa $opts{'thres_nt_switch_para'}", 0) and &stopErr("[Err] Failed to run ete3 for $tGlob{'curr_order'} -th COG\n"); 
	} else {
		&exeCmd_1cmd("$opts{'ete3_cmd'} -a ${pref}.prot.fa -o ${pref}_ete3O/", 0) and &stopErr("[Err] Failed to run ete3 for $tGlob{'curr_order'} -th COG\n"); 
	}
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
	for my $tk (qw/curr_pref curr_fh_oProt curr_toRM curr_fh_oCds/) {
		delete($glob{$tk}); 
	}
	$glob{'curr_pref'} = "$wrk_dir/$glob{'curr_order'}"; 
	&tsmsg(join('', "[Msg]   Prepared data for [$glob{'curr_pref'}] with [", scalar( map { @$_ } @{$glob{'curr_tab'}}),"] genes\n")); 
	$glob{'curr_fh_oProt'} = &openFH("$glob{'curr_pref'}.prot.fa", '>'); 
	defined $opts{'in_cds'} and $glob{'curr_fh_oCds'}  = &openFH("$glob{'curr_pref'}.cds.fa" , '>'); 
	my %h; 
	for my $id ( map { @$_ } @{$glob{'curr_tab'}} ) {
		defined $h{$id} and next; 
		defined $prot_fas{$id} or &stopErr("[Err] No prot-seq found for ID [$id]\n"); 
		print {$glob{'curr_fh_oProt'}} ">$prot_fas{$id}{'key'}\n$prot_fas{$id}{'seq'}\n"; 
		defined $opts{'in_cds'} and defined $cds_fas{$id} and print {$glob{'curr_fh_oCds'}} ">$cds_fas{$id}{'key'}\n$cds_fas{$id}{'seq'}\n"; 
	}
	delete($glob{'curr_tab'}); 
	close( $glob{'curr_fh_oProt'} ); delete( $glob{'curr_fh_oProt'} ); 
	defined $opts{'in_cds'} and do { close( $glob{'curr_fh_oCds'} ); delete( $glob{'curr_fh_oCds'} ); }; 
	push(@{$glob{'curr_toRM'}}, "$glob{'curr_pref'}.prot.fa"); 
	push(@{$glob{'file_toRM'}}, @{$glob{'curr_toRM'}}); 
	my %t_glob = %glob; 
	delete($t_glob{'file_toRM'}); 
	delete($t_glob{'sepPref'}); 
	push( @{$glob{'sepPref'}}, \%t_glob ); 
}#sub prepare_ete3()


