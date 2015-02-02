#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

use Getopt::Long; 
my %opts; 
Getoptions(\%opts, 
	"help!", 
	"faF:s", "gffF:s", 
	"scfID:s", "suff:s", 
); 
my $faF  = $opts{'faF'} // 'PG1All_v2_Scf.unmask.fa'; 
my $gffF = $opts{'gffF'} // 'r1_all.gff3'; 

my $id = $opts{'scfID'} // shift; 

my $add = $opts{'suff'} // ''; 


open F,'<',"$gffF" or die; 
open O,'>',"cur$add.gff3" or die; 
while (<F>) {
	m/^$id\t/o or next; 
	chomp; 
	my @ta = split(/\t/, $_); 
#	unless ($ta[1] =~ m/^pred_gff:augustus|maker|pred_gff:augustus_masked|pred_gff:snap_masked$/) {
#		$ta[8] =~ s!Name=[^;\s]+;?!!g; 
		$ta[8] =~ s!Target=[^;\s]+;?!!g; 
#	}
	print O join("\t", @ta)."\n"; 
}
close O; 
close F; 
exeCmd( "deal_fasta.pl $faF -res $id > cur$add.fasta" ); 
