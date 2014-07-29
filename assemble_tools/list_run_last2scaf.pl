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

for (my $i=1; $i<@list; $i++) {
	my $step_dir = "Step$i"; 
	my $step_refFa = ( $i == 1 ) ? $list[0] : $prev_scfFa ; 
	my $step_qryFa = $list[$i]; 
	my $step_oPref = "S$i"; 
	my $step_scfPref = $step_oPref; 
	my $step_dbTag = "db/dbS$i"; 
	mkdir($step_dir); 
	&tsmsg("[Rec] chdir to $step_dir\n"); 
	chdir($step_dir); 
	mkdir("db/"); 
	&tsmsg("[Rec] Copy ../$step_refFa ../$step_qryFa to $step_dir\n"); 
	system "cp -p ../$step_refFa ../$step_qryFa ./"; 
	system "bash ../run_last_to_scaffold.sh $step_refFa $step_qryFa $step_dbTag $step_oPref $step_scfPref"; 
	$prev_scfFa = "${step_oPref}.scafSeq"; 
	&tsmsg("[Rec] Copy $prev_scfFa to upper dir.\n"); 
	system "cp -p $prev_scfFa ../"; 
	&tsmsg("[Rec] chdir back to ../\n"); 
	chdir("../"); 
}

&tsmsg("[Rec] All done.\n"); 

