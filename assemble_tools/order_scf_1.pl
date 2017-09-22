#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 

# Output: 
#   5 # CG_ChrID        Final_ID        Final_Direction Ordered Annotation
#   6 # CG_Chr01        CGscf_0024      R       Y
#   7 # CG_Chr01        CGscf_0019      R       Y

my $pref = 'P1Cma_Chr'; 
@ARGV > 1 and $pref = $ARGV[1]; 
splice(@ARGV, 1, 1); 

my $min_rat = 0.8; 
my @seg; # [LG_ID, scfID, [scfPos, ...], [lgPos, ...], p_cnt, m_cnt, ttl_cnt]
my $lg_dist = 1; 
while (<>) {
	chomp; 
	m/^\s*(#|$)/ and next; 
	m/^group/ and next; 
	my @ta = split(/\t/, $_); 
	my ($scfID, $scfPos, $lgID, $lgPos, $lgID_1) = @ta; 
	if ( @seg == 0 ) {
		push(@seg, [$lgID, $scfID, [$scfPos], [$lgPos], 0, 0, 1]); 
		next; 
	}
	unless ( $seg[-1][0] eq $lgID and $seg[-1][1] eq $scfID ) {
		push(@seg, [$lgID, $scfID, [$scfPos], [$lgPos], 0, 0, 1]); 
		next; 
	}
	$lgPos < $seg[-1][3][-1] and &tsmsg("[Wrn] Be careful!!! lgPos : $lgPos < $seg[-1][3][-1] for $scfID : $seg[-1][1]\n"); 
	$seg[-1][6] ++; 
	if ( $lgPos-$seg[-1][3][-1] < $lg_dist ) {
		push(@{$seg[-1][3]}, $lgPos); 
		push(@{$seg[-1][2]}, $scfPos); 
		next; 
	}
	push(@{$seg[-1][3]}, $lgPos); 
	if ( $scfPos > $seg[-1][2][-1] ) {
		$seg[-1][4] ++; 
		push(@{$seg[-1][2]}, $scfPos); 
	} elsif ( $scfPos < $seg[-1][2][-1] ) {
		$seg[-1][5] ++; 
		push(@{$seg[-1][2]}, $scfPos); 
	} else {
		die "Repeat scfPos $scfPos\n"; 
	}
}
# my @seg; # [LG_ID, scfID, [scfPos, ...], [lgPos, ...], p_cnt, m_cnt, ttl_cnt]
print STDOUT join("\t", qw/ChrID ScfID ScfDirection LgRange/)."\n"; 
for (my $i=0; $i<@seg; $i++) {
	my @ta = @{$seg[$i]}; 
	$ta[0] =~ s!^(POP1-)?LG|^UL!!; 
	$ta[0] = sprintf("%s%02d", $pref, $ta[0]); 
	my $tmp_ttl = $ta[4]+$ta[5]; 
	my $fr = &assign_fr( $ta[4], $ta[5] ); 
	my $lg_range = "$ta[3][0]-$ta[3][-1]:T$ta[6]:f$ta[4]:m$ta[5]"; 
	print STDOUT join("\t", $ta[0], $ta[1], $fr, $lg_range)."\n"; 
}

sub assign_fr {
	my ($fN, $rN) = @_; 
	my $tN = $fN+$rN; 
	my $back = 'U'; 
	$tN <= 0 and return ($back); 
	if ( $tN != 4 ) {
		$fN >= $tN * $min_rat and $back = 'F'; 
		$rN >= $tN * $min_rat and $back = 'R'; 
	} elsif ( $tN == 4 ) {
		$fN >= 3 and $back = 'F'; 
		$rN >= 3 and $back = 'R'; 
	} else {
		; 
	}
	return($back); 
}

