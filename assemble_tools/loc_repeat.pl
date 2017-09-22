#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

!@ARGV and die "perl $0 repeat_depth_min NSP306_1kb_Cor1st_P2Gt5h.bwt2.srt.bam.dep > NSP306_1kb_Cor1st_P2Gt5h.bwt2.srt.bam.dep_repDdepth.loc_id_start_end_length\n"; 

my $repDep_min = 109; 
$repDep_min = shift; 

my $dist_rep = 100; 
# my $min_rep = 200; 

my %raw_loc; 
while (<>) {
	$. % 1000e3 == 1 and &tsmsg("[Msg] Reading $. line\n"); 
	chomp; 
	my ($scfID, $scfPos, $repNum) = split(/\t/, $_); 
	$repNum >= $repDep_min or next; 
	if ( defined $raw_loc{$scfID} ) {
		if ( $raw_loc{$scfID}[-1][1] >= $scfPos-1 ) {
			$raw_loc{$scfID}[-1][1] < $scfPos and $raw_loc{$scfID}[-1][1] = $scfPos; 
		} else {
			push( @{$raw_loc{$scfID}}, [$scfPos, $scfPos] ); 
		}
	} else {
		push(@{$raw_loc{$scfID}}, [$scfPos, $scfPos]); 
	}
}

for my $scfID ( sort keys %raw_loc ) {
	my @merged_loc; 
	for (my $i=0; $i<@{$raw_loc{$scfID}}; $i++) {
		if ($i == 0) {
			push(@merged_loc, [ @{$raw_loc{$scfID}[$i]} ]); 
			next; 
		}
		if ( $raw_loc{$scfID}[$i][0] - $merged_loc[-1][1] - 1 <= $dist_rep ) {
			$merged_loc[-1][1] < $raw_loc{$scfID}[$i][1] and $merged_loc[-1][1] = $raw_loc{$scfID}[$i][1]; 
		} else {
			push(@merged_loc, [ @{$raw_loc{$scfID}[$i]} ]); 
		}
	}
	for my $ar1 (@merged_loc) {
		print STDOUT join("\t", $scfID, $ar1->[0], $ar1->[1], $ar1->[1]-$ar1->[0]+1)."\n"; 
	}
}

