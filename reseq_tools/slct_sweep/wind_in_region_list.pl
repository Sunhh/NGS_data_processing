#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $mat_obj = mathSunhh->new(); 


my $help_txt = <<HH; 

perl $0 Spo_V1p5_w50ks5k.p2.dvd_4to1.compare.slct.merged Spo_V1p5_w50ks5k.p2.grp1_grp4.PIavg 

HH
@ARGV == 2 or &LogInforSunhh::usage($help_txt); 

my %region_list_need = %{ &load_region( $ARGV[0] ) }; 
my %region_list_in   = %{ &load_region( $ARGV[1] ) }; 

for my $chrID (sort { $region_list_in{$a}[0][3] <=> $region_list_in{$b}[0][3] } keys %region_list_in) {
	defined $region_list_need{$chrID} or next; 
	@{$region_list_in{$chrID}}   = sort { $a->[0] <=> $b->[0] } @{$region_list_in{$chrID}}; 
	@{$region_list_need{$chrID}} = sort { $a->[0] <=> $b->[0] } @{$region_list_need{$chrID}}; 
	for my $tr1 ( @{$region_list_in{$chrID}} ) {
		my $is_ok = 0; 
		for my $tr2 ( @{$region_list_need{$chrID}} ) {
			$tr2->[0] > $tr1->[1] and last; 
			$tr2->[1] < $tr1->[0] and next; 
			my @tt = $mat_obj->ovl_region( $tr1->[0], $tr1->[1], $tr2->[0], $tr2->[1] ); 
			if ( $tt[0] > 0 and $tt[0] == $tr1->[1]-$tr1->[0]+1 ) {
				$is_ok = 1; 
				last; 
			}
		}
		if ( $is_ok == 1 ) {
			print STDOUT "$tr1->[2]\n"; 
		}
	}
}


sub load_region {
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($. == 1 and $ta[0] =~ m/^(chrID|chr|chromID)$/i) {
			next; 
		}
		push(@{$back{$ta[0]}}, [$ta[1], $ta[2], $_, $.]); 
	}
	close($fh); 
	return(\%back); 
}
