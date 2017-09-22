#!/usr/bin/perl -w
use strict; 
use warnings; 
use fileSunhh; 
use LogInforSunhh; 
use mathSunhh; 
my $ms_obj = mathSunhh->new(); 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"wind_size:i", # 100 
	"wind_step:i", # wind_size
	"wind_start:s", # 1
	"wind_end:s", # Max of table. 
); 
$opts{'wind_size'} //= 100; 
$opts{'wind_step'} //= $opts{'wind_size'}; 
$opts{'wind_start'} //= ''; 
$opts{'wind_end'}   //= ''; 

my %glob; 

my $help_txt = <<HH; 
perl $0 whole_grp_2ato1/apple_whole.chr1.snp.out.xpclr.txt

-wind_size   [100]
-wind_step   [-wind_size]
-wind_start  [pos_min]
-wind_end    [pos_max]

HH

!@ARGV and &LogInforSunhh::usage($help_txt); 

my %xpclr = %{ &load_xpclr( $ARGV[0] ) }; 
$opts{'wind_start'} eq '' and $opts{'wind_start'} = $xpclr{'min'}; 
$opts{'wind_end'}   eq '' and $opts{'wind_end'}   = $xpclr{'max'}; 
my %wind = %{ $ms_obj->setup_windows('ttl_start' => $opts{'wind_start'}, 'ttl_end'=>$opts{'wind_end'}, 'wind_size'=>$opts{'wind_size'}, 'wind_step'=> $opts{'wind_step'}, 'minRatio'=> 0.8) }; 

my %score; 
my %glb; 
for my $ar1 (@{$xpclr{'line'}}) {
	$glb{'chrID'} //= $ar1->[0]; 
	my $snp_num = $ar1->[2]; 
	$snp_num > 0 or next; 
	my $p = $ar1->[3]; 
	my $v = $ar1->[5]; 
	$v =~ m/^inf$/i and next; 
	$p < $opts{'wind_start'} and next; 
	$p > $opts{'wind_end'} and next; 
	my @si_ar = @{ $ms_obj->map_windows( 'wind_hash'=>\%wind, 'posi' => $p ) }; 
	for my $si ( @si_ar ) {
		push(@{$score{$si}}, $v); 
	}
}

print STDOUT join("\t", qw/chrID WindS WindE WindLen Avg Cnt AdjAvg Min Max/)."\n"; 
for my $si (sort {$a<=>$b} keys %score) {
	my ($ws, $we, $wlen) = @{$wind{'loci'}{$si}}; 
	my %stat = %{ &mathSunhh::ins_calc( $score{$si}, 0 ) }; 
	print STDOUT join("\t", $glb{'chrID'}, $ws, $we, $wlen, @stat{qw/MEAN COUNT interval_mean min max/})."\n"; 
}


# [Sunhh@Penguin by_Map]$ head -4 whole_grp_2ato1/apple_whole.chr1.snp.out.xpclr.txt
# 1    0     3                 278.000000   0.000249    inf         0.000000
# 1    1     4                 378.000000   0.000339    inf         0.000000
# 1    2     4                 478.000000   0.000429    inf         0.000000
# chr# grid# #ofSNPs_in_window physical_pos genetic_pos XPCLR_score max_s
sub load_xpclr {
	my $fn = shift; 
	my $fh = &openFH($fn, '<'); 
	my %back; 	
	while (&wantLineC($fh)) {
		my @ta = &splitL(" ", $_); 
		$ta[3] = int($ta[3]); 
		$back{'max'} //= $ta[3]; 
		$back{'min'} //= $ta[3]; 
		$back{'max'} < $ta[3] and $back{'max'} = $ta[3]; 
		$back{'min'} > $ta[3] and $back{'min'} = $ta[3]; 
		push(@{$back{'line'}}, [@ta]); 
	}
	close($fh); 
	return(\%back); 
}
