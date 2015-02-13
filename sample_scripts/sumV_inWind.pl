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
	"colN:i", 
	"maxLine:i", 
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
#  -colN           [2]
#
#  -maxLine        [-1]
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
$opts{'colN'} = $opts{'colN'} // 2; 

## Setup windows. 
my $mm = mathSunhh->new(); 

#  CHR_A         BP_A        SNP_A  CHR_B         BP_B        SNP_B           R2
#      1         1284      s1_1284      1         1906      s1_1906    0.0801282
my $inLine = 0; 
my %chr_wind; 
while (<>) {
	chomp; m/^\s*$/ and next; 
	my @ta = split(/\s+/, $_); 
	$ta[0] eq 'chr' and next; 
	defined $chr_wind{$ta[0]} or $chr_wind{$ta[0]} = $mm->setup_windows(
	  'ttl_start'=>$opts{'wind_start'}, 
	  'ttl_end' => $opts{'wind_end'}, 
	  'wind_size' => $opts{'wind_length'}, 
	  'wind_step' => $opts{'wind_step'}, 
	  'minRatio' => 0 
	); 
	
	$opts{'maxLine'} > 0 and $inLine > $opts{'maxLine'} and last; 
	my $dist = $ta[1]; 
	my (@wind_i) = @{ $mm->map_windows( 'posi' => $dist , 'wind_hash' => $chr_wind{$ta[0]} ) }; 
	for my $ti ( @wind_i ) {
		push(@{$chr_wind{$ta[0]}{'vv'}{$ti}}, $ta[$opts{'colN'}]); 
	}
	$inLine ++; 
}

print STDOUT join("\t", qw/ChromID WindS WindE SUM COUNT MEAN/)."\n"; 
for my $chrID ( sort keys %chr_wind ) {
	for my $ti ( @{$chr_wind{$chrID}{'info'}{'windSloci'}} ) {
		defined $chr_wind{$chrID}{'vv'}{$ti} or next; 
		my %ins_back = %{ mathSunhh::ins_calc( $chr_wind{$chrID}{'vv'}{$ti} ) }; 
		for my $tk (qw/MEAN/) {
			$ins_back{$tk} = sprintf("%.4f", $ins_back{$tk}); 
		}
		print STDOUT join("\t", $chrID, @{$chr_wind{$chrID}{'loci'}{$ti}}[0,1], @ins_back{qw/SUM COUNT MEAN/})."\n"; 
	}
}


