#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"wind_len:i", # default 100
	"wind_step:i", # Default -wind_len 
); 

my $help_txt = <<HH; 

perl $0 pen_R01.ctg.fa.chop_info in.bam.dep > pen_R01.ctg.fa.chop_info_aD

# Method to generate in.chop_info : 
# perl deal_fasta.pl pen_R01.scf_Gt5h.fa -chop_seq -chop_len 100 -chop_step 100 -chop_min 100 -chop_info > pen_R01.scf_Gt5h.fa.chop_info
# ### The .chop_info format : key     WI      WS      WE      GC      AG      Wkey    All     N    wN

-help

-wind_len       [100]. 
-wind_step      [-wind_len]. 


HH

$opts{'wind_len'} //= 100; 
$opts{'wind_step'} //= $opts{'wind_len'}; 

!@ARGV and &LogInforSunhh::usage($help_txt); 

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
	my ( $idx_aref ) = &idx_by_WL_WStep( $ta[1], $opts{'wind_len'}, $opts{'wind_step'} ); 
	for my $idx ( @$idx_aref ) {
		defined $wInfo{$ta[0]}{$idx} or next; 
		push(@{$wInfo{$ta[0]}{$idx}[1]}, $ta[2]); 
	}
}

sub idx_by_WL_WStep {
	my ($p, $len, $step, $start) = @_; 
	$start //= 1; 

	my @back; 
	my $med_idx = int( ($p-$start+1-1)/$len )+1; 
	for (my $i=$med_idx; $i>0; $i--) {
		$start+($i-1)*$step+$len-1 < $p and last; 
		$start+($i-1)*$step > $p and &stopErr("[Err] Something wrong! $start+($i-1)*$step > $p (\$p)\n"); 
		unshift(@back, $i); 
	}
	for (my $i=$med_idx+1; $start+($i-1)*$step <= $p; $i++) {
		push(@back, $i); 
	}

	return(\@back); 
} # sub idx_by_WL_WStep() 

for my $cid (sort keys %wInfo) {
	for my $wi ( sort { $a <=> $b } keys %{$wInfo{$cid}} ) {
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

