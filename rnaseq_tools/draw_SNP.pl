#!/usr/bin/perl -w 
use strict; 
use SVG; 
use Getopt::Long; 
use LogInforSunhh; 
my %opts; 
GetOptions(\%opts, 
	"bpPerPoint:f", # 30000 ; 
	"lis_chrLen:s", # ITAG2.3_genomic.fa.kl , required. 
	"lis_windR:s",  # CheesRIL_p204_LA483_diffHomo_wi35offs.tab.indv_cntR 
	"indv_ID:s",    # individual ID. 
	"P1_ID:s",      # Parent 1 ID (blue color)
	"P2_ID:s",      # Parent 2 ID (red color)
	"out_svg:s",    # out.svg 

	"min_siteN:i",  # 0 
); 

# Format of ITAG2.3_genomic.fa.kl : 
#  key     len
#  SL2.40ch00      21805821
#  SL2.40ch01      90304244
#  SL2.40ch02      49918294

# Format of CheesRIL_p204_LA483_diffHomo_wi35offs.tab.indv_cntR : 
#  chr             binID   wind_S  wind_E  P1_R    P2_R    Hete_R  Miss_R  sumN
#  SL2.40ch00      2       400001  600000  0.00    0.00    100.00  0.00    1
#  SL2.40ch00      5       1000001 1200000 42.86   14.29   42.86   0.00    7
#  SL2.40ch00      28      5600001 5800000 100.00  0.00    0.00    0.00    4

my $help_txt = <<HH; 

perl $0 
  -lis_chrLen      ITAG2.3_genomic.fa.kl 
  -lis_windR       CheesRIL_p204_LA483_diffHomo_wi35offs.tab.indv_cntR 
  -out_svg         out.svg 

  -indv_ID         individual ID
  -P1_ID           P1
  -P2_ID           P2

  -bpPerPoint      30000
  -min_siteN       0

HH

( defined $opts{'lis_chrLen'} and defined $opts{'lis_windR'} ) or &LogInforSunhh::usage($help_txt); 

$opts{'bpPerPoint'} //= 30e3; 
$opts{'min_siteN'}  //= 0; 
$opts{'out_svg'}    //= 'out.svg'; 

my %chrL = %{ &load_chrLen($opts{'lis_chrLen'}) }; 
my %cntR = %{ &load_wind_cntR($opts{'lis_windR'}) }; 
my $indvID = ( defined $opts{'indv_ID'} ) ? $opts{'indv_ID'} : $opts{'lis_windR'} ; 
my $P1_ID = ( defined $opts{'P1_ID'} ) ? $opts{'P1_ID'} : 'P1'; 
my $P2_ID = ( defined $opts{'P2_ID'} ) ? $opts{'P2_ID'} : 'P2'; 

my $chr_longest = (sort {$b <=> $a} values %chrL)[0]; 
my $chr_number = scalar( keys %chrL ); 

my %plot_par; 
$plot_par{'chr_height'}   //= 40; 
$plot_par{'chr_distance'} //= 20; 

$plot_par{'width'}  //= $chr_longest / $opts{'bpPerPoint'} + 50; 
$plot_par{'height'} //= 50 + ($plot_par{'chr_height'}+$plot_par{'chr_distance'}+4) * $chr_number; 
$plot_par{'y_0'}    //= $plot_par{'height'} - 50; 
$plot_par{'x_0'}    //= 200; 

$plot_par{'colWid'} //= 0.5; 
$plot_par{'font_size'} //= 15; 

$plot_par{'y_height_percent'} = $plot_par{'chr_height'} / 100; 

$plot_par{'grp_color'} = {
'P1_R' =>  'blue', 
'P2_R' =>  'red',
'Hete_R' => 'green', 
'Miss_R' => 'grey'
}; 

my $svg = SVG->new( 'width'=>$plot_par{'x_0'} + $plot_par{'width'}, 'height'=>$plot_par{'height'} ); 
my %grps; 
$grps{'P1_R'} = $svg->group (
	'id' => 'P1_R', 
	'stroke' => $plot_par{'grp_color'}{'P1_R'}, 'fill' => $plot_par{'grp_color'}{'P1_R'}, 
	'stroke-width' => $plot_par{'colWid'}, 
); 

$grps{'P2_R'} = $svg->group (
	'id' => 'P2_R', 
	'stroke' => $plot_par{'grp_color'}{'P2_R'}, 'fill' => $plot_par{'grp_color'}{'P2_R'}, 
	'stroke-width' => $plot_par{'colWid'}, 
); 

$grps{'Hete_R'} = $svg->group (
	'id' => 'Hete_R', 
	'stroke' => $plot_par{'grp_color'}{'Hete_R'}, 'fill' => $plot_par{'grp_color'}{'Hete_R'}, 
	'stroke-width' => $plot_par{'colWid'}, 
); 

$grps{'Miss_R'} = $svg->group (
	'id' => 'Miss_R', 
	'stroke' => $plot_par{'grp_color'}{'Miss_R'}, 'fill' => $plot_par{'grp_color'}{'Miss_R'}, 
	'stroke-width' => $plot_par{'colWid'}, 
); 
$grps{'chr_bone'} = $svg->group (
	'id' => 'chr_bone', 
	'stroke' => 'black', 'fill' => 'none', 
	'stroke-width' => 1, 
	'text-anchor' => 'middle', 
	'font-weight' => 'normal', 
	'font-size'   => $plot_par{'font_size'}, 
	'font-family' => 'Arial', 
); 

$grps{'chr_bone'}->text(
  'x' => $plot_par{'x_0'}, 'y' => 20, 
  -cdata => "$indvID", 
  'text-anchor' => 'start', 
); 
{
	my $cum_x = 0; 
	my %tr2tb = qw( Hete_R Heterozygous Miss_R Missing ); 
	$tr2tb{'P1_R'} = $P1_ID; 
	$tr2tb{'P2_R'} = $P2_ID; 
	for my $tr (qw/P1_R P2_R Hete_R Miss_R/) {
		defined $tr2tb{ $tr } or die "unknown tr [$tr]\n"; 
		my $tb = $tr2tb{ $tr }; 
		$grps{$tr}->rectangle(
		  'x'  => $plot_par{'x_0'}+$cum_x, 'y' => 40, 
		  'width' => 10, 'height' => 10, 
		); 
		$grps{$tr}->text(
		  'x'  => $plot_par{'x_0'}+$cum_x+13, 'y' => 50, 
		  -cdata => "$tb", 
		  'text-anchor' => 'start', 
		); 
		$cum_x += (10+13+80); 
	}
}

$plot_par{'cumu_y'} = $plot_par{'y_0'}; 
for my $chrn (reverse sort keys %chrL) {
	# Chromosome box : 
	$grps{'chr_bone'}->rectangle(
	  'x' => $plot_par{'x_0'}, 'y' => $plot_par{'cumu_y'} - $plot_par{'chr_height'}, 
	  'width' => $chrL{$chrn}/$opts{'bpPerPoint'}, 'height' => $plot_par{'chr_height'}
	); 

	$grps{'chr_bone'}->text(
	  'x' => $plot_par{'x_0'}-20, 'y' => $plot_par{'cumu_y'} - $plot_par{'chr_height'}/2, 
	  -cdata => "CHR: $chrn", 
	  'text-anchor' => 'end', 
	); 

	for (my $i=0; $i<=$chrL{$chrn}; $i+=5e6) {
		my $j = sprintf("%0.1f", $i/1e6); 
		$grps{'chr_bone'}->text(
		  'x' => $plot_par{'x_0'} + $i/$opts{'bpPerPoint'}, 'y' => $plot_par{'cumu_y'}+10, 
		  -cdata => "$j M", 
		  'text-anchor' => 'middle', 
		  'font-weight' => 'normal', 
		  'font-size'   => 10, 
		); 
	}

	for my $ta ( @{$cntR{$chrn}} ) {
		my %vR; 
		@vR{qw/ws we P1_R P2_R Hete_R Miss_R sumN/} = @$ta; 
		my $cum_delt = 0; 
		for my $tb (qw/P1_R P2_R Hete_R Miss_R/) {
			$vR{$tb} > 0 or next; 
			my $th = $vR{$tb} * $plot_par{'y_height_percent'}; 
			$grps{$tb}->rectangle(
			  'x' => $plot_par{'x_0'} + $vR{'ws'}/$opts{'bpPerPoint'}, 'y'=>$plot_par{'cumu_y'}-$cum_delt-$th, 
			  'width' => ($vR{'we'}-$vR{'ws'}+1)/$opts{'bpPerPoint'}, 'height' => $th, 
			); 
			$cum_delt += $th; 
		}
	}
	$plot_par{'cumu_y'} -= ($plot_par{'chr_height'}+$plot_par{'chr_distance'}); 
}
open OS,'>',"$opts{'out_svg'}" or die; 
print OS $svg->xmlify; 
close OS; 




###### Sub-routines. 

sub load_wind_cntR {
	my $fn = shift; 
	my %back; 
	open F,'<',"$fn" or die; 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[0] eq 'chr' and next; 
		$ta[8] >= $opts{'min_siteN'} or next; 
		push(@{$back{$ta[0]}}, [@ta[2,3, 4,5,6,7, 8]]); 
	}
	close F; 
	return(\%back); 
}

sub load_chrLen {
	my $fn = shift; 
	my %back; 
	open F,'<',$fn or die; 
	while (<F>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		$ta[0] eq 'key' and next; 
		$back{$ta[0]} = $ta[1]; 
	}
	close F; 
	return (\%back); 
}


