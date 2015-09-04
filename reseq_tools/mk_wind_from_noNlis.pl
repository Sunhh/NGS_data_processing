#!/usr/bin/perl
use strict; 
use warnings; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 

@ARGV == 3 or &usage(); 

sub usage {
print <<HH;

perl $0 wind_length wind_step Malus_x_domestica.v1.0-primary.pseudo.fa.noN_list

# 2015-06-02 Count non-N bp number within windows given by [wind_length wind_step]. 
#  Format of Malus_x_domestica.v1.0-primary.pseudo.fa.noN_list: 
#   Key     Length  MatchStart      MatchEnd        MatchLen
#   chr10   38388735        1       5022            5022
#   chr10   38388735        5024    10768           5745
HH
exit 0; 
}

my $wind_len = shift ; 
my $wind_step = shift ; 

my %wind; 
my %w2v; 
my %si2idx; 
my %idx2si; 
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($chrID, $chrLen, $mS, $mE, $mL) = @ta; 
	$chrLen =~ m/^length$/i and next; 
	defined $wind{$chrID} or $wind{$chrID} = $ms_obj->setup_windows( 'ttl_start'=>1, 'ttl_end'=>$chrLen, 'wind_size'=>$wind_len, 'wind_step'=>$wind_step ); 
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
		$w2v{$chrID}{$si} += ($ovlE-$ovlS+1); 
	}
}

print STDOUT join("\t", qw/ChromID WindS WindE WindL BpCnt/)."\n"; 
for my $chrID ( sort keys %wind ) {
	for my $si ( sort { $a<=>$b } keys %{$wind{$chrID}{'loci'}} ) {
		my $vv = $w2v{$chrID}{$si} // 0; 
		my ($wS, $wE, $iL) = @{$wind{$chrID}{'loci'}{$si}}; 
		print STDOUT join("\t", $chrID, $wS, $wE, $iL, $vv)."\n"; 
	}
}


