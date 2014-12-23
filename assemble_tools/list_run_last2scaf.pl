#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

my @list; 
while (<>) {
	chomp; 
	m/^\s*$/ and next; 
	push(@list, (split(/\s+/, $_))[0]); 
}
&tsmsg("[Rec]total ", $#list+1," files to be joined.\n"); 

my $prev_scfFa = ''; 

my %glob_tools; 
$glob_tools{cp} = "cp"; 

for (my $i=1; $i<@list; $i++) {
	my $step_dir = "Step$i"; 
	my $step_refFa = ( $i == 1 ) ? $list[0] : $prev_scfFa ; 
	( defined $step_refFa and $step_refFa ne '' and -f "$step_refFa" ) or &stopErr("[Err] Please check refFa=[$step_refFa]\n"); 
	my $step_qryFa = $list[$i]; 
	my $step_oPref = "S${i}"; 
	my $step_scfPref = "${step_oPref}_"; 
	my $step_dbTag = "db/dbS$i"; 
	mkdir($step_dir); 
	&tsmsg("[Rec] chdir to $step_dir\n"); 
	chdir($step_dir); 
	mkdir("db/"); 
	&exe_cmd( "$glob_tools{cp} -p ../$step_refFa ../$step_qryFa ./" ); 
	$prev_scfFa = &run_last2scaff( $step_refFa, $step_qryFa, $step_dbTag, $step_oPref, $step_scfPref ); 
	&exe_cmd("$glob_tools{cp} -p $prev_scfFa ../"); 
	&tsmsg("[Rec] chdir back to ../\n"); 
	chdir("../"); 
}

&tsmsg("[Rec] All done.\n"); 

sub run_last2scaff {
	my ( $refFa, $qryFa, $dbTag, $oPref, $scfPref ) = @_; 
	( defined $refFa and -f "$refFa" ) or &stopErr("[Err] Please check refFa=[$refFa]\n"); 
	( defined $qryFa and -f "$qryFa" ) or &stopErr("[Err] Please check qryFa=[$qryFa]\n"); 
	defined $dbTag or $dbTag = "dbTag"; 
	defined $oPref or $oPref = "oPref"; 
	defined $scfPref or $scfPref = "scfPref"; 

	# Medium file names. 
	my %med; 
	$med{refClip} = "${refFa}.clip"; 
	$med{qryClip} = "${qryFa}.clip"; 
	$med{refCtg}  = "${oPref}_Ref.fitMugsy.fa"; 
	$med{qryCtg}  = "${oPref}_Qry.fitMugsy.fa"; 
	$med{ctgFa}   = "${oPref}_all.fitMugsy.fa"; 
	$med{rawMaf}  = "${oPref}.maf"; 
	$med{fitMaf}  = "${oPref}.fitMugsy.maf"; 
	$med{xmfa}    = "${oPref}.fitMugsy.maf.xmfa"; 
	$med{pref_mugO} = "${oPref}_syn"; 
	$med{link}      = "$med{pref_mugO}.paired.maf.link"; 
	$med{tbl}       = "${oPref}.tbl"; 
	$med{scf}       = "${oPref}.scf"; 
	$med{scfSeq}    = "${oPref}.scfSeq"; 
	$med{drop}      = "$med{scfSeq}.dropped"; 

	# Exe tools: 
	my %tools; 
	my $HOME = `echo \$HOME`; 
	chomp($HOME); 
	my $SRC  = "/data/Sunhh/src"; 

	$tools{pl_clipScaf}    = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/clip_scaf_end.pl"; 
	$tools{pl_addTag2fa}   = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/add_tag_to_fsa.pl"; 
	$tools{pl_formatMAF}   = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/format_maf_forMugsy.pl"; 
	$tools{pl_maf2fasta}   = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/maf2fasta.pl"; 
	$tools{pl_pairedMAF}   = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/get_paired_maf.pl"; 
	$tools{pl_goodLink}    = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/good_link_fromMaf.pl"; 
	$tools{pl_joinScaf}    = "perl $HOME/tools/github/NGS_data_processing/assemble_tools/join_link.pl"; 
	$tools{exe_lastdb}     = "$HOME/bin/lastdb"; 
	$tools{exe_lastal}     = "$HOME/bin/lastal"; 
	$tools{exe_paraFa}     = "$HOME/bin/parallel-fasta"; 
	$tools{exe_lastSp}     = "$HOME/bin/last-split"; 
	$tools{sh_mugsyenv}    = "source $SRC/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyenv.sh"; 
	$tools{exe_mugsyWGA}   =        "$SRC/Align/Mugsy/mugsy_x86-64-v1r2.3/mugsyWGA"; 
	$tools{cat}            = "cat"; 

	# Carry out commands. 
	&exe_cmd("$tools{pl_clipScaf} $refFa > $med{refClip}"); 
	&exe_cmd("$tools{pl_clipScaf} -min_outLen 500 $qryFa > $med{qryClip}"); 

	&exe_cmd("$tools{pl_addTag2fa} Ref $med{refClip} > $med{refCtg}"); 
	&exe_cmd("$tools{pl_addTag2fa} Qry $med{qryClip} > $med{qryCtg}"); 
	&exe_cmd("$tools{cat} $med{refCtg} $med{qryCtg} > $med{ctgFa}"); 

	&exe_cmd("$tools{exe_lastdb} -w 5 -m111111111111111111111111111110 $dbTag $med{refCtg}"); 
	&exe_cmd("$tools{exe_paraFa} \"$tools{exe_lastal} -q3 -e80 $dbTag | $tools{exe_lastSp} \" < $med{qryCtg} > $med{rawMaf}"); 
	&exe_cmd("$tools{pl_formatMAF} -specs \",\" $med{rawMaf} -out $med{fitMaf}"); 
	&exe_cmd("$tools{pl_maf2fasta} $med{fitMaf} > $med{xmfa}"); 

	&exe_cmd("$tools{sh_mugsyenv} ; $tools{exe_mugsyWGA} --outfile $med{pref_mugO} --seq $med{ctgFa} --aln $med{xmfa} --distance 2000 --minlength 50 1>stdout.$med{pref_mugO} 2>stderr.$med{pref_mugO}"); 

	&exe_cmd("$tools{pl_pairedMAF} $med{pref_mugO}.maf > $med{pref_mugO}.paired.maf"); 
	&exe_cmd("$tools{pl_goodLink} $med{pref_mugO}.paired.maf -minClustLen 2000 1>$med{link} 2>$med{link}.err"); 
	&exe_cmd("$tools{pl_joinScaf} $med{link} -inCtgFa $med{ctgFa} -outLnkTbl $med{tbl} -outScfLnk $med{scf} -outScfFas $med{scfSeq} -scfPref $scfPref -dropScfIn $med{drop}"); 

	&tsmsg("[Rec] Finish last2scaff [ $refFa, $qryFa, $dbTag, $oPref, $scfPref ]\n"); 

	return $med{scfSeq}; 
}#End sub run_last2scaff


sub exe_cmd {
	for (@_) {
		chomp($_); 
		&tsmsg("[CMD] $_\n"); 
		system("$_"); 
	}
	return 0; 
}







