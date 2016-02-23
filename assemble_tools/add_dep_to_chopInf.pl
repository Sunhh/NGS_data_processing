#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 

!@ARGV and die "perl $0 pen_R01.ctg.fa.chop_info in.bam.dep > pen_R01.ctg.fa.chop_info_aD\n"; 

my $info_file = shift; 

my %wInfo = %{ &load_chop_info($info_file) }; 

# [Sunhh@panda 01.from_JC]$ head -4 P201512_PG1_toCtg_primary.noClip_nmR0p01_refGt3h.bam.dep
# 6954706 11      1
# 6954706 12      1
# 6954706 13      1
while (<>) {
	$. % 10e6 == 1 and &tsmsg("[Msg] $. line.\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	defined $wInfo{$ta[0]} or next; 
	my $idx = int( ($ta[1]-1)/100 )+1; 
	defined $wInfo{$ta[0]}{$idx} or next; 
	push(@{$wInfo{$ta[0]}{$idx}[1]}, $ta[2]); 
}

for my $cid (keys %wInfo) {
	for my $wi ( keys %{$wInfo{$cid}} ) {
		if ( @{$wInfo{$cid}{$wi}[1]} < 50 ) {
			print STDOUT "$wInfo{$cid}{$wi}[0]\t0\t0\n"; 
			next; 
		}
		my %cc = %{ &mathSunhh::ins_calc( $wInfo{$cid}{$wi}[1], 1 ) }; 
		my $avg1 = sprintf("%0.4f", $cc{'SUM'} / 100 ); 
		my $avg2 = sprintf("%0.4f", $cc{'interval_mean'} * $cc{'COUNT'} / 100 ); 
		print STDOUT "$wInfo{$cid}{$wi}[0]\t$avg1\t$avg2\n"; 
	}
}

sub load_chop_info {
	my ($fn) = @_; 
	# key             WI      WS      WE      GC      AG      Wkey            len
	# 136463376       1       1       100     0.28    0.74    136463376_1     100
	# 136463378       1       1       100     0.34    0.75    136463378_1     100
	# 136463380       1       1       100     0.19    0.99    136463380_1     100
	# 136463382       1       1       100     0.3     0.64    136463382_1     100
	my %back; 
	my $fh = &openFH($fn,'<'); 
	my %cnt; 
	{
		my $hd = <$fh>; 
		chomp($hd); 
		print STDOUT "$hd\tAvgDep1\tAvgDep2\n"; # 
	}
	while (<$fh>) {
		$cnt{'line'} ++; 
		$cnt{'line'} % 1e6 == 1 and &tsmsg("[Msg] $cnt{'line'} in $fn\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		my ($key, $wi) = @ta[0,1]; 
		$back{$key}{$wi} = [$_, [], $cnt{'line'}]; 
	}
	close($fh); 
	return(\%back); 
}

