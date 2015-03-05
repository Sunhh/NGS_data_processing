#!/usr/bin/perl
use strict; 
use warnings; 
use LogInforSunhh; 
use mathSunhh; 

use Getopt::Long; 
my %opts;
GetOptions(\%opts, 
	"help!", 
	"wind_length:i", "wind_start:i", "wind_end:i", "wind_step:i", 
	"chr_colN:i", "pos_colN:i", "cnt_colN:i", 
	"maxLine:i", "showAll!", "trimTail!", 
	"skipHN:i", 
); 

-t and !@ARGV and &usage(); 
$opts{'help'} and &usage(); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 input_plink.ld
#
#  -help
#
#  -wind_length    [1]
#  -wind_step      [wind_length]
#  -wind_start     [1]
#  -wind_end       [9999999]
#
#  -chr_colN       [0] 999999 means there is no chr_colN being used. 
#  -pos_colN       [1] 999999 means there is no pos_colN being used, then setup windows from cnt_colN
#  -cnt_colN       [] Count will be always 1 if not given this column number. 
#
#  -showAll        [] Show all windows if given. 
#  -trimTail       [] Trim tailing windows. 
#
#  -maxLine        [-1]
#  -skipHN         [0] Skip header lines number. 
#
#  -symbol         ['\t']
################################################################################
HH
	exit 1; 
}

# chr10   8955    0.198064517483858       0.15073145245559        -0.281317715822552      CT      T       T       CT      T       T       T       T       T
# chr10   9225    0.19718189638096        0.468664169787765       1.6131216677027 CT      T       T       T       T       T       C       C       T       T
# chr10   9318    0.196326356382151       0.217152412804587       0.123718584278465       C       C       C       C       CT      C       AC      C       C
#

$opts{'wind_length'} = $opts{'wind_length'} // 1; 
$opts{'wind_start'} = $opts{'wind_start'} // 1; 
$opts{'wind_end'} = $opts{'wind_end'} // 9999999; 
$opts{'min_r2'} = $opts{'min_r2'} // 0; 
$opts{'wind_step'} = $opts{'wind_step'} // $opts{'wind_length'}; 
$opts{'maxLine'} = $opts{'maxLine'} // -1; 
$opts{'chr_colN'} = $opts{'chr_colN'} // 0 ;
$opts{'pos_colN'} = $opts{'pos_colN'} // 1 ;
$opts{'skipHN'} = $opts{'skipHN'} // 0; 
my $symbol = "\t"; 
defined $opts{'symbol'} and $symbol = $opts{'symbol'}; 
# $opts{'cnt_colN'} = $opts{'cnt_colN'} ;

## Setup windows. 
my $mm = mathSunhh->new(); 

#  CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2
#      1         1284      s1_1284      1         1906      s1_1906    0.0801282
&tsmsg("[Rec] Reading file.\n"); 
my $inLine = 0; 
my %chr_wind; 
my $ln_cnt = 0; 
while (<>) {
	$ln_cnt ++; 
	$ln_cnt > $opts{'skipHN'} or next; 
	$. % 1e6 == 1 and &tsmsg("[Msg] Reading [$.] line(s).\n"); 
	chomp; m/^\s*$/ and next; 
	my @ta = split(/$symbol/o, $_); 
	my $chrV = ( $opts{'chr_colN'} == 999999 ) ? 0     : $ta[ $opts{'chr_colN'} ] ; 
	$chrV =~ m/^chr$/i and next; 
	my $cntV = ( defined $opts{'cnt_colN'}  ) ? $ta[ $opts{'cnt_colN'} ] : 1 ; 
	my $posV = ( $opts{'pos_colN'} == 999999 ) ? $cntV : $ta[ $opts{'pos_colN'} ] ; 
	( defined $cntV and $cntV ne '' ) or next; 
	defined $chr_wind{$chrV} or $chr_wind{$chrV} = $mm->setup_windows(
	  'ttl_start' => $opts{'wind_start'}, 
	  'ttl_end'   => $opts{'wind_end'}, 
	  'wind_size' => $opts{'wind_length'}, 
	  'wind_step' => $opts{'wind_step'}, 
	  'minRatio'  => 0 
	); 
	
	$opts{'maxLine'} > 0 and $inLine > $opts{'maxLine'} and last; 
	my (@wind_i) = @{ $mm->map_windows( 'posi' => $posV , 'wind_hash' => $chr_wind{$chrV} ) }; 
	for my $ti ( @wind_i ) {
		push(@{$chr_wind{$chrV}{'vv'}{$ti}}, $cntV); 
	}
	$inLine ++; 
}
&tsmsg("[Rec] Writing table.\n"); 

my %endIdx; 
if ( $opts{'trimTail'} ) {
	for my $chrID ( sort keys %chr_wind ) {
		my @idx = @{$chr_wind{$chrID}{'info'}{'windSloci'}} ; 
		my $last_i = $#idx; 
		for ( ; $last_i >= 0 ; $last_i -- ) {
			defined $chr_wind{$chrID}{'vv'}{ $idx[$last_i] } and last; 
		}
		$endIdx{$chrID} = ($last_i == -1) ? 'stop' : $idx[$last_i] ; 
	}
}

print STDOUT join("\t", qw/ChromID WindS WindE SUM COUNT MEAN/)."\n"; 
for my $chrID ( sort keys %chr_wind ) {
	for my $ti ( @{$chr_wind{$chrID}{'info'}{'windSloci'}} ) {
		$opts{'showAll'} or defined $chr_wind{$chrID}{'vv'}{$ti} or next; 
		$opts{'trimTail'} and ( $endIdx{$chrID} eq 'stop' or $endIdx{$chrID} < $ti ) and last; 
		defined $chr_wind{$chrID}{'vv'}{$ti} or $chr_wind{$chrID}{'vv'}{$ti} = []; 
		my %ins_back = %{ mathSunhh::ins_calc( $chr_wind{$chrID}{'vv'}{$ti} ) }; 
		for my $tk (qw/MEAN/) {
			$ins_back{$tk} ne '' and $ins_back{$tk} = sprintf("%.4f", $ins_back{$tk}); 
		}
		print STDOUT join("\t", $chrID, @{$chr_wind{$chrID}{'loci'}{$ti}}[0,1], @ins_back{qw/SUM COUNT MEAN/})."\n"; 
	}
}
&tsmsg("[Rec] All done.\n"); 


