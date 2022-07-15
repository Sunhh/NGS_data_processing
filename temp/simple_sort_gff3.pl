#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 WM97_v1.chr.annot_prot.gff3 > WM97_v1.chr.annot_prot.gff3.srt\n"; 

my %l1_gene; 
my %gene_to_mRNA; 
my %l2_mRNA; 
my %l3_CDS; 

while (<>) {
	m!^\s*#|^\s*$! and next; 
	chomp; 
	my @ta=split(/\t/, $_); 
	@ta = @ta; 
	my $id = ''; 
	$ta[8] =~ m!\bID=([^\s;]+)! and $id = $1; 
	my $pid = ''; 
	$ta[8] =~ m!\bParent=([^\s;]+)! and $pid = $1; 
	if ($ta[2] eq 'gene') {
		$id eq '' and die "$_\n"; 
		defined $l1_gene{$id} and die "repG: $_\n"; 
		$l1_gene{$id} = [@ta]; 
	} elsif ($ta[2] eq 'mRNA' or $ta[2] eq 'transcript') {
		$id eq '' and die "tid:$_\n"; 
		$pid eq '' and die "tidP:$_\n"; 
		defined $l2_mRNA{$id} and die "repM: $_\n"; 
		$l2_mRNA{$id} = [@ta]; 
		push(@{$gene_to_mRNA{$pid}}, $id); 
	} elsif ($ta[2] eq 'CDS' or $ta[2] eq 'exon') {
		$pid eq '' and die "cidP:$_\n"; 
		push(@{$l3_CDS{$pid}}, [@ta]); 
	} else {
		die "$_\n"; 
	}
}
for my $gid (sort {$l1_gene{$a}[0] cmp $l1_gene{$b}[0] || $l1_gene{$a}[3] <=> $l1_gene{$b}[3] || $l1_gene{$a}[4] <=> $l1_gene{$a}[4]} keys %l1_gene) {
	print STDOUT join("\t", @{$l1_gene{$gid}})."\n"; 
	defined $gene_to_mRNA{$gid} or next; 
	for my $mid (sort { $l2_mRNA{$a}[3] <=> $l2_mRNA{$b}[3] || $l2_mRNA{$a}[4] <=> $l2_mRNA{$b}[4] } @{$gene_to_mRNA{$gid}}) {
		print STDOUT join("\t", @{$l2_mRNA{$mid}})."\n"; 
		defined $l3_CDS{$mid} or next; 
		if ($l2_mRNA{$mid}[6] eq '+') {
			@{$l3_CDS{$mid}} = sort { $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] } @{$l3_CDS{$mid}}; 
		} elsif ($l2_mRNA{$mid}[6] eq '-') {
			@{$l3_CDS{$mid}} = sort { $b->[4] <=> $a->[4] || $b->[3] <=> $a->[4] } @{$l3_CDS{$mid}}; 
		} else {
			die "@{$l2_mRNA{$mid}}\n"; 
		}
		for my $cline (@{$l3_CDS{$mid}}) {
			print STDOUT join("\t", @{$cline})."\n"; 
		}
	}
}


