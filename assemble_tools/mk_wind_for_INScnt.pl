#!/usr/bin/perl
use strict; 
use warnings; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

@ARGV == 4 or &usage(); 

sub usage {
print <<HH;

perl $0 INS_avg_Len wind_length wind_step Malus_x_domestica.v1.0-primary.pseudo.fa.noN_list

# 2015-06-02 Count non-N bp number within windows given by [wind_length wind_step]. 
#  Format of Malus_x_domestica.v1.0-primary.pseudo.fa.noN_list: 
#   Key     Length  MatchStart      MatchEnd        MatchLen
#   chr10   38388735        1       5022            5022
#   chr10   38388735        5024    10768           5745
HH
exit 0; 
}

my $ins_len = shift ; 
my $wind_len = shift ; 
my $wind_step = shift ; 

my %wind; 
my %w2v; 
my %si2idx; 
my %idx2si; 
my %si_to_mSE; # si -to- mSE_region_in_window_in_noN_list $si_to_mSE{$chrID}{$start_position} => [ [$mS, $mE], [$mS, $mE], ... ]; 
my %chr_len; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($chrID, $chrLen, $mS, $mE, $mL) = @ta; 
	$chrLen =~ m/^length$/i and next; 
	defined $wind{$chrID} or $wind{$chrID} = $ms_obj->setup_windows( 'ttl_start'=>1, 'ttl_end'=>$chrLen, 'wind_size'=>$wind_len, 'wind_step'=>$wind_step ); 
	$chr_len{$chrID} //= $chrLen; 
	unless ( defined $si2idx{$chrID} ) {
		my $idx = 0; 
		for my $si ( sort { $a<=>$b } keys %{$wind{$chrID}{'loci'}} ) {
			$si2idx{$chrID}{$si} = $idx; 
			$idx2si{$chrID}{$idx} = $si; 
			$idx ++; 
		}
	}
	my @si_mS = sort { $a <=> $b } @{ $ms_obj->map_windows( 'posi'=>$mS, 'wind_hash'=>$wind{$chrID} ) }; 
	my @si_mE = sort { $a <=> $b } @{ $ms_obj->map_windows( 'posi'=>$mE, 'wind_hash'=>$wind{$chrID} ) }; 
	for my $idx ( $si2idx{$chrID}{ $si_mS[0] } .. $si2idx{$chrID}{ $si_mE[-1] } ) {
		my $si = $idx2si{$chrID}{$idx}; 
		my ($wS, $wE, $iL) = @{ $wind{$chrID}{'loci'}{$si} }; 
		my $ovlS = &mathSunhh::max($wS, $mS); 
		my $ovlE = &mathSunhh::min($wE, $mE); 
		push( @{$si_to_mSE{$chrID}{$si}}, [$ovlS, $ovlE] ); 
		$w2v{$chrID}{$si}[0] += ($ovlE-$ovlS+1); 
	}
}

for my $chrID (sort keys %si_to_mSE) {
	for my $si ( sort { $a <=> $b } keys %{$si_to_mSE{$chrID}} ) {
		my (@fwd_blks, @rev_blks); 
		for my $from_se (@{ $si_to_mSE{$chrID}{$si} }) {
			my ($from_s, $from_e) = @$from_se; 
			my ($fwd_s, $fwd_e) = ( $from_s+$ins_len-1, $from_e+$ins_len-1 ); 
			$fwd_e > $chr_len{$chrID} and $fwd_e = $chr_len{$chrID}; 
			if ( $fwd_s > $chr_len{$chrID} ) {
				; 
			} else {
				push( @fwd_blks, @{ &blks_by_region( $fwd_s, $fwd_e, $chrID ) } ); 
			}
			my ($rev_s, $rev_e) = ( $from_s-$ins_len+1, $from_e-$ins_len+1 ); 
			$rev_s > 0 or $rev_s = 1; 
			if ( $rev_e < 1 ) {
				; 
			} else {
				push( @rev_blks, @{ &blks_by_region( $rev_s, $rev_e, $chrID ) } ); 
			}
		}
		my (@merged_fwd_blk, @merged_rev_blk); 
		scalar(@fwd_blks) > 0 and @merged_fwd_blk = @{ $ms_obj->mergeLocBlk( \@fwd_blks ) }; 
		scalar(@rev_blks) > 0 and @merged_rev_blk = @{ $ms_obj->mergeLocBlk( \@rev_blks ) }; 
		$w2v{$chrID}{$si}[1] //= 0; # For good forward length 
		$w2v{$chrID}{$si}[2] //= 0; # For good reverse length 
		for my $se (@merged_fwd_blk) {
			$w2v{$chrID}{$si}[1] += ($se->[1]-$se->[0]+1); 
		}
		for my $se (@merged_rev_blk) {
			$w2v{$chrID}{$si}[2] += ($se->[1]-$se->[0]+1); 
		}
	}
}

print STDOUT join("\t", qw/ChromID WindS WindE WindL BpCnt FwdBp RevBp/)."\n"; 
for my $chrID ( sort keys %wind ) {
	for my $si ( sort { $a<=>$b } keys %{$wind{$chrID}{'loci'}} ) {
		$w2v{$chrID}{$si}[0] //= 0; 
		$w2v{$chrID}{$si}[1] //= 0; 
		$w2v{$chrID}{$si}[2] //= 0; 
		my ($wS, $wE, $iL) = @{$wind{$chrID}{'loci'}{$si}}; 
		print STDOUT join("\t", $chrID, $wS, $wE, $iL, $w2v{$chrID}{$si}[0], $w2v{$chrID}{$si}[1], $w2v{$chrID}{$si}[2])."\n"; 
	}
}


sub blks_by_region {
	my ($mS, $mE, $chrID) = @_; 
	my @back_blks; 
	my @si_mS = sort { $a<=>$b } @{ $ms_obj->map_windows( 'posi'=>$mS, 'wind_hash'=>$wind{$chrID} ) }; 
	my @si_mE = sort { $a<=>$b } @{ $ms_obj->map_windows( 'posi'=>$mE, 'wind_hash'=>$wind{$chrID} ) }; 
	for my $idx ( $si2idx{$chrID}{ $si_mS[0] } .. $si2idx{$chrID}{ $si_mE[-1] } ) {
		my $si = $idx2si{$chrID}{$idx}; 
		defined $si_to_mSE{$chrID}{$si} or next; 
		for my $ref_se (@{ $si_to_mSE{$chrID}{$si} }) {
			my ($ref_s, $ref_e) = @$ref_se; 
			$ref_s > $mE and next; 
			$ref_e < $mS and next; 
			my $ovlS = &mathSunhh::max($ref_s, $mS); 
			my $ovlE = &mathSunhh::min($ref_e, $mE); 
			push(@back_blks, [$ovlS, $ovlE]); 
		}
	}
	return \@back_blks; 
}

