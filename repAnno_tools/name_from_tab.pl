#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

-t and !@ARGV and die "perl $0 dgt.tab\n"; 

# eleID   eleS    eleE    Str     seqID   LTR1_S  LTR1_E  LTR2_S  LTR2_E  Inner_S Inner_E PBS_S   PBS_E   PPT_S   PPT_E   scfID
# RR1     98334   105001  ?       seq10   98339   99198   104137  104996  99199   104136  -1      -1      -1      -1      S400016_pilon
# RR2     785204  788754  +       seq10   785209  785369  788589  788749  785370  788588  785373  785384  -1      -1      S400016_pilon

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ( $ta[0] eq 'eleID' ) {
		# print STDOUT "$_\n"; 
		next; 
	}
	my ($inner_s, $inner_e, $pbs_s, $pbs_e, $ppt_s, $ppt_e) = @ta[9,10, 11,12, 13,14]; 
	my ($eleID, $seqID, $scfID) = @ta[0,4,15]; 
	my $tk1 = "${eleID}_${seqID}_$ta[1]_$ta[2]_ELE_$scfID"; # Element region. (with TSD)
	my $tk2 = "${eleID}_${seqID}_$ta[5]_$ta[8]_LTR_$scfID"; # LTR region. (without TSD)
	my $tk3 = "${eleID}_${seqID}_$ta[9]_$ta[10]_INN_$scfID"; # Internal region. (without ltr region)
	print STDOUT "$tk1\t${eleID}_\t${eleID}\n"; 
	print STDOUT "$tk2\t${eleID}_\t${eleID}\n"; 
	print STDOUT "$tk3\t${eleID}_\t${eleID}\n"; 
}

