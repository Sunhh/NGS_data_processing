#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 

-t and !@ARGV and &stopErr("perl $0 aug1.out > aug1_withID_for_gff3_merge.gff\n"); 

my %is_good; 
%is_good = qw(
	gene         gene
	transcript   mRNA
	mrna         mRNA
	cds          CDS
); 

my @geneLines; 
my $mID = ''; 
while (<>) {
	m/^\s*$/ and next; 
	m/^\s*#/ and next; 
	m/^>/ and last; 
	chomp; 
	my @ta = split(/\t/, $_); 
	$ta[2] = lc($ta[2]); 
	defined $is_good{$ta[2]} or next; 
	$ta[2] = $is_good{$ta[2]}; 
	$ta[2] eq 'gene' and next; 
	if ( $ta[2] eq 'mRNA' ) {
		if ( scalar(@geneLines) > 0 ) {
			&fmt_gL(\@geneLines); 
			&out_gL(\@geneLines); 
			@geneLines = (); 
		}
		$ta[8] =~ m/^(g\d+\.t\d+)$/ or &stopErr("[Err] 2: ta8=$ta[8]\n"); 
		$mID = $1; 
		push(@geneLines, [ [@ta], $mID ]); 
	} elsif ( $ta[2] eq 'CDS' ) {
		$ta[8] =~ m/^transcript_id \"([^"]+)\"; gene_id \"([^"]+)\";?$/ or &stopErr("[Err] 3: ta8=$ta[8]\n"); 
		my $t_mID = $1; 
		$mID eq '' and $mID = $t_mID; 
		$mID eq $t_mID or &stopErr("[Err] 4: Changed mID [$mID] to [$t_mID]\n"); 
		push(@geneLines, [ [@ta], $t_mID ]); 
	} else {
		&stopErr("[Err] Why here\n"); 
	}
}

if ( scalar(@geneLines) > 0 ) {
	&fmt_gL(\@geneLines); 
	&out_gL(\@geneLines); 
	@geneLines = (); 
}

sub fmt_gL {
	my $ar = shift; 
	my $mID = ''; 
	for my $r1 (@$ar) {
		$mID eq '' or last; 
		$r1->[1] ne '' and $mID = $r1->[1]; 
	}
	my $cdsN = 0; 
	for my $r1 (@$ar) {
		if ( $r1->[0][2] eq 'mRNA' ) {
			$r1->[0][2] = 'match'; 
			$r1->[0][8] = "ID=$mID;"; 
		} elsif ( $r1->[0][2] eq 'CDS' ) {
			$r1->[0][2] = 'match_part'; 
			$cdsN ++; 
			$r1->[0][8] = "ID=${mID}.c$cdsN;Parent=$mID;"; 
		} else {
			&stopErr("[Err] Why here.\n"); 
		}
	}
	return 0; 
}

sub out_gL {
	my $ar = shift; 
	for my $r1 (@$ar) {
		print STDOUT join("\t", @{$r1->[0]})."\n"; 
	}
	return 0; 
}

