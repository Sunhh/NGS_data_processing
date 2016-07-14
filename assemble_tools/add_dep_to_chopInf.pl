#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"wind_len:i", # default 100
	"wind_step:i", # Default -wind_len 
	"slct_chr:s@", # Default [], Format: chrID:start-end; 
); 

my $help_txt = <<HH; 

perl $0 pen_R01.ctg.fa.chop_info in.bam.dep > pen_R01.ctg.fa.chop_info_aD

# Method to generate in.chop_info : 
# perl deal_fasta.pl pen_R01.scf_Gt5h.fa -chop_seq -chop_len 100 -chop_step 100 -chop_min 100 -chop_info > pen_R01.scf_Gt5h.fa.chop_info
# ### The .chop_info format : key     WI      WS      WE      GC      AG      Wkey    All     N    wN

-help

-wind_len       [100]. 
-wind_step      [-wind_len]. 

-slct_chr       [] chrID:start-end ....


HH

$opts{'wind_len'} //= 100; 
$opts{'wind_step'} //= $opts{'wind_len'}; 

!@ARGV and &LogInforSunhh::usage($help_txt); 
$opts{'help'} and &LogInforSunhh::usage($help_txt); 

my %chr_windows; 

my %slct_chr; 
my $has_slct_chr = 0; 
if ( defined $opts{'slct_chr'} and scalar(@{$opts{'slct_chr'}}) > 0 ) {
	$has_slct_chr = 1; 
	for my $a1 ( @{$opts{'slct_chr'}} ) {
		$a1 =~ s!\s!!g; 
		my ($id, $s, $e); 
		if ( $a1 =~ m/^(\S+):(\d+)\-(\d+)$/ ) {
			($id, $s, $e) = ($1, $2, $3); 
		} elsif ( $a1 =~ m/^(\S+)$/ ) {
			$id = $1; 
			$s = $e = 0; 
		} else {
			&stopErr("[Err] Bad format for -slct_chr : [$a1]\n"); 
		}
		push( @{$slct_chr{$id}}, [$s, $e] ); 
	}
}

my $info_file = shift; 

my %wInfo = %{ &load_chop_info($info_file) }; 
my %wMM; 
for my $cid (keys %wInfo) {
	for my $wi ( keys %{$wInfo{$cid}} ) {
		my @ta = split(/\t/, $wInfo{$cid}{$wi}[0]); 
		$wMM{$cid}{'max'} //= $ta[3]; 
		$wMM{$cid}{'min'} //= $ta[2]; 
		$wMM{$cid}{'max'} < $ta[3] and $wMM{$cid}{'max'} = $ta[3]; 
		$wMM{$cid}{'min'} > $ta[2] and $wMM{$cid}{'min'} = $ta[2]; 
	}
}

# [Sunhh@panda 01.from_JC]$ head -4 P201512_PG1_toCtg_primary.noClip_nmR0p01_refGt3h.bam.dep
# 6954706 11      1
# 6954706 12      1
# 6954706 13      1
my %tc = ( 'cntN_base'=>0 , 'cntN_step'=>5e6 ); 
while (<>) {
	&fileSunhh::log_section( $., \%tc ) and &tsmsg("[Msg] $. in .dep file\n"); 
	chomp; 
	my @ta = split(/\t/, $_); 
	defined $wInfo{$ta[0]} or next; 
	( $ta[1] >= $wMM{$ta[0]}{'min'} and $ta[1] <= $wMM{$ta[0]}{'max'} ) or next; 
	my ( $idx_aref ) = &idx_by_WL_WStep( $ta[1], $opts{'wind_len'}, $opts{'wind_step'} ); 
	for my $idx ( @$idx_aref ) {
		defined $wInfo{$ta[0]}{$idx} or next; 
		push(@{$wInfo{$ta[0]}{$idx}[1]}, $ta[2]); 
	}
}

for my $cid (sort keys %wInfo) {
	for my $wi ( sort { $wInfo{$cid}{$a}[2] <=> $wInfo{$cid}{$b}[2] } keys %{$wInfo{$cid}} ) {
		if ( @{$wInfo{$cid}{$wi}[1]} < 1 ) {
			# There is no reads mapping here. 
			print STDOUT "$wInfo{$cid}{$wi}[0]\t0\t0\t0\t0\n"; 
			next; 
		}
		my %cc = %{ &mathSunhh::ins_calc( $wInfo{$cid}{$wi}[1], 1 ) }; 
		my $avg1 = ( $wInfo{$cid}{$wi}[3] > 0 ) ? ($cc{'SUM'} / $wInfo{$cid}{$wi}[3]) : 0 ;
		$avg1 = sprintf("%0.4f", $avg1 ); 
		my $avg2 = sprintf("%0.4f", $cc{'interval_mean'} ); # Instead of averaged by all good sites, I use covered sites here. 
		my $avg3 = sprintf("%0.4f", $cc{'MEAN'}); 
		print STDOUT "$wInfo{$cid}{$wi}[0]\t$avg1\t$avg2\t$cc{'SUM'}\t$avg3\n"; 
	}
}

sub load_chop_info {
	my ($fn) = @_; 
	# key             WI      WS      WE      GC      AG      Wkey            len
	# 136463376       1       1       100     0.28    0.74    136463376_1     100
	# 136463378       1       1       100     0.34    0.75    136463378_1     100
	# 136463380       1       1       100     0.19    0.99    136463380_1     100
	# 136463382       1       1       100     0.3     0.64    136463382_1     100
	# key             WI      WS      WE      GC      AG      Wkey            len             N       wN
	# Cma_Chr01       1       1       10000   0.4181  0.4787  Cma_Chr01_1     13080099        149064  852
	# Cma_Chr01       2       5001    15000   0.3981  0.4738  Cma_Chr01_2     13080099        149064  852
	# Cma_Chr01       3       10001   20000   0.3881  0.4715  Cma_Chr01_3     13080099        149064  0
	#
	my %back; 
	my $fh = &openFH($fn,'<'); 
	my %cnt; 
	{
		my $hd = <$fh>; 
		chomp($hd); 
		print STDOUT "$hd\tAvgDep\tAdjAvgCovDep\tSumDep\tAvgCovDep\n"; # 
	}
	my %tc = ( 'cntN_base'=>0 , 'cntN_step'=>1e6 ); 
	while (<$fh>) {
		$cnt{'line'} ++; 
		&fileSunhh::log_section( $cnt{'line'}, \%tc ) and &tsmsg("[Msg] $cnt{'line'} in $fn\n"); 
		chomp; 
		my @ta = split(/\t/, $_); 
		my ($key, $wi) = @ta[0,1]; 
		my $wlen = $ta[3]-$ta[2]+1; 
		if ( defined $ta[9] ) {
			$wlen = $wlen - $ta[9]; 
		}

		if ( $has_slct_chr == 1 ) {
			defined $slct_chr{ $key } or next; 
			my $is_in = 0; 
			for my $a1 (@{$slct_chr{$key}}) {
				( $a1->[0] == 0 or $a1->[0] <= $ta[3] ) or next; 
				( $a1->[1] == 0 or $a1->[1] >= $ta[2] ) or next; 
				$is_in = 1; 
				last; 
			}
			$is_in == 1 or next; 
		}

		unless ( defined $chr_windows{$key} ) {
			($chr_windows{$key}) = $ms_obj->setup_windows(
			  'ttl_start' => 1, 
			  'ttl_end'   => $ta[7], # Chr length 
			  'wind_size' => $opts{'wind_len'}, 
			  'wind_step' => $opts{'wind_step'}, 
			); 
		}

		if ( defined $chr_windows{$key}{'info'}{'windSloci'}[$wi-1] ) {
			my $si = $chr_windows{$key}{'info'}{'windSloci'}[$wi-1]; 
			my $ei = $chr_windows{$key}{'loci'}{$si}[1]; 
			( $si == $ta[2] and $ei == $ta[3] ) or &stopErr("[Err] -wind_len $opts{'wind_len'} -wind_step $opts{'wind_step'} doesn't agree with data line : $_\n"); 
		}

		$back{$key}{$wi} = [$_, [], $cnt{'line'}, $wlen]; 
	}
	close($fh); 
	return(\%back); 
}

sub idx_by_Pos {
	my ($p, $wind_hr) = @_; 
	my @back; 

	my $si_ar = $ms_obj->map_windows( 'position' => $p, 'wind_hash' => $wind_hr ); 
	for my $si (@$si_ar) {
		push(@back, $wind_hr->{'info'}{'Sloc2wi'}{$si}); 
	}

	return(\@back); 
}

sub idx_by_WL_WStep {
	my ($p, $len, $step, $start) = @_; 
	$start //= 1; 

	my @back; 
	my $med_idx = int( ($p-$start+1-1)/$step )+1; 
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

