#!/usr/bin/perl
# Only written for maker_gff3. 
use LogInforSunhh; 
use strict; 
use warnings; 

-t and !@ARGV and &stopErr("perl $0 mRNA_name.list maker.gff3\n"); 

my $lisF = shift; 
my $gffF = shift; 

open F1,'<',"$lisF" or die; 
my %useM; 
my %useG; 
# maker-S401991_pilon-est_gff_Cufflinks-gene-0.0-mRNA-1
while (<F1>) {
	s/\s+$//; 
	my @ta = split(/\t/, $_); 
	my $gName = $ta[0]; 
	$gName =~ m/\s/ and &stopErr("[Err] spaced gene name [$gName]\n"); 
	$useM{$gName} = 1; 
	$gName =~ s/\-mRNA\-\d+//; 
	$useG{$gName} = 1; 
}
close F1; 
open F2,'<',"$gffF" or die; 
my $is_fa = 0; 
while (<F2>) {
	s/[^\S\t]+$//; 
	m/^\s*$/ and next; 
	my $is_out = 0; 
	if ( $is_fa == 1) {
		$is_out = 1; 
	} elsif ( m/^\s*\>/ ) {
		$is_fa = 1; 
		$is_out = 1; 
	} elsif ( m/^#/ ) {
		$is_out = 1; 
	} else {
		my @ta = split(/\t/, $_); 
		$ta[8] =~ m/^ID=([^\s;]+);/ or &stopErr("[Err] Unknown ta8=[$ta[8]]\n"); 
		my $gName = $1; 
		$gName =~ m/,/ and &stopErr("[Err] Unparsed ta8=[$ta[8]]\n"); 
		$gName =~ s/:\S*$//; 
		if ($ta[2] eq 'gene') {
			defined $useG{$gName} and $is_out = 1; 
		} else {
			defined $useM{$gName} and $is_out = 1; 
		}
	}
	$is_out == 1 and print "$_\n"; 
}
close F2; 
