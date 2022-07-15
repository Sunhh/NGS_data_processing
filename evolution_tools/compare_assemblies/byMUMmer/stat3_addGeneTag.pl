#!/usr/bin/perl
use strict; 
use warnings; 

# ==> WM97_v1.scf.annot_prot.gff3 <==
# WM97_scaffold3460	GLEAN	gene	102	557	0.919626	+	.	ID=Cla000001
# WM97_scaffold3460	GLEAN	mRNA	102	557	0.919626	+	.	ID=Cla000001.t1;Parent=Cla000001;
# WM97_scaffold3460	GLEAN	CDS	102	557	.	+	0	Parent=Cla000001.t1
# WM97_scaffold8822	GLEAN	gene	65	295	0.831796	+	.	ID=Cla000002
# 
# ==> ngs.nlis.coverTag <==
# Key	Length	MatchStart	MatchEnd	MatchLen	Covered
# WM97_scaffold11	169061	3288	3336	49	1
# WM97_scaffold11	169061	14316	14502	187	1
# WM97_scaffold11	169061	58876	58908	33	1

-t and !@ARGV and die "perl $0 WM97_v1.scf.annot_prot.gff3 ngs.nlis.coverTag > ngs.nlis.coverTag.geneTag\n"; 

my $gffF = shift; 
my $nlis = shift; 

my $dist = 500; 

open F1,'<',"$gffF" or die; 
my %loc; 
while (<F1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] eq 'mRNA' or next; 
	if ( $ta[6] eq '+' ) {
		my $min = $ta[3]-$dist; 
		$min > 0 or $min = 1; 
		push(@{$loc{$ta[0]}}, [$min, $ta[4]]); 
	} elsif ( $ta[6] eq '-' ) {
		push(@{$loc{$ta[0]}}, [$ta[3], $ta[4]+$dist]); 
	}
}
close F1; 
for my $k1 (keys %loc) {
	@{$loc{$k1}} = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]  } @{$loc{$k1}}; 
}

open F2,'<',"$nlis" or die; 
while (<F2>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[0] eq 'Key'){
		print join("\t", @ta, "GeneTag")."\n"; 
		next; 
	}
	my $is_gene = 0; 
	for my $k1 (@{$loc{$ta[0]}}) {
		$k1->[1] < $ta[2] and next; 
		$k1->[0] > $ta[3] and last; 
		$is_gene = 1; 
		last; 
	}
	print join("\t", @ta, $is_gene)."\n"; 
}
close F2; 

