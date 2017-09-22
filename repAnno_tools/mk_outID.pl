#!/usr/bin/perl
use strict; 
use warnings; 
# seqID   LTR1_S  LTR1_E  LTR2_S  LTR2_E  Inner_S Inner_E PBS_S   PBS_E   Strand  scfID
# seq10   785209  785369  788589  788749  785370  788588  785373  785384  +       S400016_pilon

while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	print STDOUT join("\t", $_, "$ta[0]_$ta[1]_$ta[4]_$ta[10]")."\n"; 
}
