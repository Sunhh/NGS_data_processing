#!/usr/bin/perl -w
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 in_prot.fa out_prot.fa\n"; 

my $fin = shift; 
my $fout = shift; 

&exeCmd_1cmd("cat $fin | deal_fasta.pl -rmTailX_prot | deal_fasta.pl -frag 0-0 -frag_width 100 -frag_head | deal_fasta.pl -chopKey ':1\\-\\d+' | deal_fasta.pl -rmDefinition > $fout"); 



