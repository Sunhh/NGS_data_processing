#!/usr/bin/perl
use strict; 
use warnings; 
use SVG; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
use mcsSunhh; 
use plotSunhh; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	# Output file : 
	"out_svg:s", # finshed. 
	# Input alignments : finished. 
	"in_tab:s", # result from "follow_mcscan.pl -aln2table". 
	# Input gene positions : MCScan formatted .gff . finished. 
	"in_pos:s", # Format : Chr_ID  \\t  Gene_ID  \\t  Gene_Start  \\t  Gene_End (unstranded)
	# Input X-Y axis : finished. 
	"chrLis_x:s", "chrLis_y:s", # Format : Chr_ID \\t Chr_Len \\n
	# For frame options : finished. 
	"frame_setopt:s@", # 
	# For syntenic blocks : finished. 
	"blk_setopt:s@", # Can be assigned multiple times. 
	# For shadow lines : finished. 
	"xAxis_shadow_setopt:s@", "yAxis_shadow_setopt:s@", 
	# For additional lines : finished. 
	"add_lineX:s", "add_lineY:s", # Add border lines for x/y axises. format : Chr_ID \\t Chr_Position \\n ; The border line will be drawn at Chr_ID - Chr_Position ; 
	# Other settings : finished. 
	"img_setopt:s@", # Can be assigned multiple times. $default_string{'img_setopt'}
	"other_setopt:s@", # Can be anything I want to setup . 
	# Timing 
	"show_time!", 
); 

################################################################################
#    Setup basic parameters. 
################################################################################
$opts{'show_time'} and &tsmsg("[Msg] Script [$0] begin.\n"); 

my $outSvgFh = \*STDOUT ; 
defined $opts{'out_svg'} and $outSvgFh = &openFH( $opts{'out_svg'}, '>' ); 


my %default_string; 
$default_string{'img_setopt'} = ['title=;horizMargin=50;vertMargin=50;xAxisLabel=Genome 1;yAxisLabel=Genome 2']; 
$default_string{'blk_setopt'} = [
	'l_lwd=1;l_col=black;',                                                # lines linking points within block 
	'lks_col_min=0;lks_col_max=3;lks_col_rgb=blue:cyan:green:yellow:red;', # block's average Ks : ''
	'p_radius=1;p_lwd=0;p_col=black',                                      # points standing for gene pairs : ''
	'pks_col_min=0;pks_col_max=3;pks_col_rgb=blue:cyan:green:yellow:red;', # gene pairs' Ks     : ''
]; 
$default_string{'frame_setopt'} = [
	'title_FontSize=40;title_FontFam=ArialNarrow;title_FontWeight=bold;title_HoriAln=middle;title_VerAln=text-after-edge;', 
	
	'scaleBar_height=10;scaleBar_width=80;scaleBar_x0=5;scaleBar_y0=5;scaleBar_dvdN=100', 
	
	'xAxis_TickStep=10000000;xAxis_TickUnit=M;xAxis_TickLwd=1;xAxis_TickLen=3;xAxis_TickCol=blue;xAxis_TickFontSize=0;xAxis_FontWeight=lighter;xAxis_FontFam=ArialNarrow;xAxis_TxtHoriAln=middle;xAxis_TxtVertAln=text-before-edge;', 
	# Axis_TickLwd   <= 0 means there is no ticks to plot. 
	# Axis_TickUnit could be 'bp/K/M/none'. 'bp/K/M' is for bp/Kb/Mb, and 'none' ignores the plot of ticks' names. 
	'xAxis_BpPoint=300000;xAxis_BorderLwd=1;xAxis_BorderCol=black;xAxis_SepBorder=0;', 
	'xAxis_LabFontSize=20;xAxis_LabFontFam=ArialNarrow;xAxis_LabFontWeight=bold;xAxis_LabHoriAln=middle;xAxis_LabVertAln=text-before-edge;', 
	'xAxis_ChrIDFontSize=5', 
	# Axis_BorderLwd <= 0 means there is no border to plot. 
	# If Axis_SepBorder > 0, there will be an Axis_SepBorder length interval between two adjacent chromosomes' borders. 
	'yAxis_TickStep=10000000;yAxis_TickUnit=M;yAxis_TickLwd=1;yAxis_TickLen=3;yAxis_TickCol=blue;yAxis_TickFontSize=0;yAxis_FontWeight=lighter;yAxis_FontFam=ArialNarrow;yAxis_TxtHoriAln=middle;yAxis_TxtVertAln=text-before-edge;', 
	'yAxis_BpPoint=300000;yAxis_BorderLwd=1;yAxis_BorderCol=black;yAxis_SepBorder=0;', 
	'yAxis_LabFontSize=20;yAxis_LabFontFam=ArialNarrow;yAxis_LabFontWeight=bold;yAxis_LabHoriAln=middle;yAxis_LabVertAln=text-before-edge;', 
	'yAxis_ChrIDFontSize=5', 
]; 
$default_string{'other_setopt'} = [
	'shadow_stroke=black;shadow_strokeWid=0.1;shadow_fill=none;', 
]; 

my (%img_opt); &_img_setopt(); 
my (%blk_opt); &_blk_setopt(); 
my (%frame_opt); &_frame_setopt(); 
my (%shad_opt); &_shadow_setopt(); # Will read in data at the same time. 
# Example 'shadow_01=file_name_1:black:10;recentDup_MaMa:file_name_2:red:10;'
#   shadow_01 file format : ChromeID \\t ChromeStart \\t ChromeEnd \\n
#   'black|red' is the filling color of rectangle for shadows; 
#   '10' is the height of horizontal rectangle (X-axis) or the width of vertical rectangle (Y-axis). 
#   The border-stroke of these rectangle is zero. 
# Structure of %shad_opt : {'xAxis_shadow_setopt|yAxis_shadow_setopt'} => [ [shadow_01, file_name_1, black, 10, { chrID_1=>[[S,E], [S,E], ...] , chrID_2=>[[S,E], [S,E], ...] }], [recentDup_MaMa, file_name_2, red, 10, {}], ... ];
my %other_opt; &_other_setopt(); # 

$opts{'help'} and &usage(); 
defined $opts{'in_tab'} or &usage(); 

sub usage {

	my %txt; 
	for my $k (keys %default_string) {
		$txt{$k} = join("\n", map { "#  $_" } @{$default_string{$k}})."\n#  "; 
	}
	print <<HH; 
################################################################################
#  perl $0 
#   -help
#
#  Required : 
#   -in_tab      out.collinearity.tab 
#   -in_pos      mcscanX.gff . Format : Chr_ID  \\t  Gene_ID  \\t  Gene_Start  \\t  Gene_End (unstranded)
#   -chrLis_x    X-Axis chromosomes. Format : Chr_ID \\t Chr_Len \\n
#   -chrLis_y    Y-Axis chromosomes. Format : Chr_ID \\t Chr_Len \\n
#
#  Optional : 
#   -out_svg         Output SVG to this file instead of STDOUT . 
#   -img_setopt      Parameters for the whole image 
$txt{'img_setopt'}
#   -frame_setopt    Parameters for frame . 
$txt{'frame_setopt'}
#   -blk_setopt      Parameters for syntenic blocks .  
$txt{'blk_setopt'}
#   -xAxis_shadow_setopt  Parameters for shadows along X-Axis . 
#   -yAxis_shadow_setopt  Parameters for shadows along Y-Axis . 
#                    Example 'shadow_01=file_name_1:black:10;recentDup_MaMa:file_name_2:red:10;'
#                      shadow_01 file format : ChromeID \\t ChromeStart \\t ChromeEnd \\n
#                      'black|red' is the filling color of rectangle for shadows; 
#                      '10' is the height of horizontal rectangle (X-axis) or the width of vertical rectangle (Y-axis). 
#                      The border-stroke of these rectangle is zero. 
# 
#   -add_lineX       The border line will be drawn at Chr_ID - Chr_Position
#   -add_lineY
#                    File format : Chr_ID \\t Chr_Position \\n
# 
#   -other_setopt    Other parameters. 
$txt{'other_setopt'}
#
#   -show_time       [Boolean]
################################################################################
HH
	exit(1); 
}

################################################################################
#    Read in data files. 
################################################################################
my @tabInfo = @{ &mcsSunhh::_readInAlnTbl( $opts{'in_tab'} ) }; 
my @tabHeader; $tabInfo[0][0] =~ m/^BlkID$/i and do { @tabHeader=@{$tabInfo[0]}; shift(@tabInfo); }; 
my ( %chr2gen, %gen2loc ); 
{ my ($a, $b) = &mcsSunhh::_readInMcsGff( $opts{'in_pos'} ); %chr2gen = %$a; %gen2loc = %$b; }
my %chrLisX = %{ &mcsSunhh::_readInChrLis( $opts{'chrLis_x'} ) }; 
my %chrLisY = %{ &mcsSunhh::_readInChrLis( $opts{'chrLis_y'} ) }; 
my (%addLineX, %addLineY); 
defined $opts{'add_lineX'} and %addLineX = %{ &_readInAddLine( $opts{'add_lineX'} ) }; 
defined $opts{'add_lineY'} and %addLineY = %{ &_readInAddLine( $opts{'add_lineY'} ) }; 


################################################################################
#    Plot frame part 
################################################################################
##############################
# SVG whole image : 
##############################
my $frame_width = &mathSunhh::max( 
	&mcsSunhh::chrP_to_plotP( 
		'chrLis' => \%chrLisX, 'chrID'=>$chrLisX{'arr'}[-1][0], 'chrP'=>$chrLisX{'arr'}[-1][1], 
		'beginPlotP' => 0, 
		'BpPoint' => $frame_opt{'xAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'xAxis_SepBorder'}, 
	) 
); 
my $frame_height = &mathSunhh::max( 
	&mcsSunhh::chrP_to_plotP( 
		'chrLis' => \%chrLisY, 'chrID'=>$chrLisY{'arr'}[-1][0], 'chrP'=>$chrLisY{'arr'}[-1][1], 
		'beginPlotP' => 0, 
		'BpPoint' => $frame_opt{'yAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'yAxis_SepBorder'}, 
	)
); 

my $svg_width  = $img_opt{'horizMargin'} * 2 + $frame_width ; 
defined $img_opt{'svg_width'}  and $svg_width = $img_opt{'svg_width'}; 
my $svg_height = $img_opt{'vertMargin'} * 2 + $frame_height ; 
defined $img_opt{'svg_height'} and $svg_height = $img_opt{'svg_height'}; 

my $svg = SVG->new( 'width' => $svg_width, 'height' => $svg_height ); 

##############################
# Define SVG groups - 
##############################
my %grps; 
$grps{'title_txt'} = $svg->group(
	'id'                => "title_txt", 
	'stroke-width'      => 0, 
	'stroke'            => 'none', 
	'opacity'           => 1, 
	'text-anchor'       => $frame_opt{'title_HoriAln'}, 
	'font-weight'       => $frame_opt{'title_FontWeight'}, 
	'font-size'         => $frame_opt{'title_FontSize'}, 
	'font-family'       => $frame_opt{'title_FontFam'}, 
	'alignment-baseline'=> $frame_opt{'title_VerAln'}, 
); 
if ( $frame_opt{'xAxis_TickLwd'} > 0 ) {
	$grps{'xAxis_Tick'} = $svg->group(
		'id'                => "xAxis_Tick", 
		'stroke-width'      => $frame_opt{'xAxis_TickLwd'}, 
		'stroke'            => $frame_opt{'xAxis_TickCol'}, 
		'opacity'           => 1, 
		'text-anchor'       => $frame_opt{'xAxis_TxtHoriAln'}, 
		'font-weight'       => $frame_opt{'xAxis_FontWeight'}, 
		'font-size'         => $frame_opt{'xAxis_TickFontSize'}, 
		'font-family'       => $frame_opt{'xAxis_FontFam'}, 
		'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
	); 
}
if ( $frame_opt{'yAxis_TickLwd'} > 0 ) {
	$grps{'yAxis_Tick'} = $svg->group(
		'id'                => "yAxis_Tick", 
		'stroke-width'      => $frame_opt{'yAxis_TickLwd'}, 
		'stroke'            => $frame_opt{'yAxis_TickCol'}, 
		'opacity'           => 1, 
		'text-anchor'       => $frame_opt{'yAxis_TxtHoriAln'}, 
		'font-weight'       => $frame_opt{'yAxis_FontWeight'}, 
		'font-size'         => $frame_opt{'yAxis_TickFontSize'}, 
		'font-family'       => $frame_opt{'yAxis_FontFam'}, 
		'alignment-baseline'=> $frame_opt{'yAxis_TxtVertAln'}, 
	); 
}
$grps{'xAxis_Label'} = $svg->group(
	'id'                => "xAxis_Label", 
	'stroke-width'      => 0, 
	'stroke'            => 'none', 
	'opacity'           => 1, 
	'text-anchor'       => $frame_opt{'xAxis_LabHoriAln'}, 
	'font-weight'       => $frame_opt{'xAxis_LabFontWeight'}, 
	'font-size'         => $frame_opt{'xAxis_LabFontSize'}, 
	'font-family'       => $frame_opt{'xAxis_LabFontFam'}, 
	'alignment-baseline'=> $frame_opt{'xAxis_LabVertAln'}, 
); 
$grps{'yAxis_Label'} = $svg->group(
	'id'                => "yAxis_Label", 
	'stroke-width'      => 0, 
	'stroke'            => 'none', 
	'opacity'           => 1, 
	'text-anchor'       => $frame_opt{'yAxis_LabHoriAln'}, 
	'font-weight'       => $frame_opt{'yAxis_LabFontWeight'}, 
	'font-size'         => $frame_opt{'yAxis_LabFontSize'}, 
	'font-family'       => $frame_opt{'yAxis_LabFontFam'}, 
	'alignment-baseline'=> $frame_opt{'yAxis_LabVertAln'}, 
); 
if ( $frame_opt{'xAxis_BorderLwd'} > 0 ) {
	$grps{'xAxis_Border'} = $svg->group(
		'id'                => "xAxis_Border", 
		'stroke-width'      => $frame_opt{'xAxis_BorderLwd'}, 
		'stroke'            => $frame_opt{'xAxis_BorderCol'}, 
		'opacity'           => 1, 
	); 
}
if ( $frame_opt{'yAxis_BorderLwd'} > 0 ) {
	$grps{'yAxis_Border'} = $svg->group(
		'id'                => "yAxis_Border", 
		'stroke-width'      => $frame_opt{'yAxis_BorderLwd'}, 
		'stroke'            => $frame_opt{'yAxis_BorderCol'}, 
		'opacity'           => 1, 
	); 
}
$grps{'blk_line'} = $svg->group(
	'id'                => "blk_line", 
	'stroke-width'      => $blk_opt{'l_lwd'}, 
	'stroke'            => $blk_opt{'l_col'}, 
	'opacity'           => 1, 
	'stroke-linecap'    => 'round', 
); 
$grps{'blk_point'} = $svg->group(
	'id'                => "blk_point", 
	'stroke-width'      => $blk_opt{'p_lwd'}, 
	'stroke'            => 'black', 
	'opacity'           => 1, 
); 

##############################
# s01 - Main frame and color scale bars 
##############################
my $frame_lx = $img_opt{'horizMargin'}; # left-X 
my $frame_by = $svg_height - $img_opt{'vertMargin'}; # bottom-Y
defined $img_opt{'frame_lx'} and $frame_lx = $img_opt{'frame_lx'}; 
defined $img_opt{'frame_by'} and $frame_by = $img_opt{'frame_by'}; 
my $frame_rx = $frame_lx + $frame_width  ; # right-X
my $frame_ty = $frame_by - $frame_height ; # top-Y
$svg->rectangle(
	'x'      => $frame_lx, 'y' => $frame_ty, 
	'width'  => $frame_width, 
	'height' => $frame_height, 
	'id'     => "Main_Frame", 
	'stroke' => 'black', 
	'fill'   => 'none', 
); 

if (@{$blk_opt{'lks_col_rgb'}} > 1 and $blk_opt{'lks_col_min'} < $blk_opt{'lks_col_max'}) {
	my $stepX = $frame_opt{'scaleBar_width'} / $frame_opt{'scaleBar_dvdN'}; 
	for (my $i=0; $i<=$frame_opt{'scaleBar_dvdN'}; $i++) {
		my $pX = $frame_opt{'scaleBar_x0'} + $frame_opt{'scaleBar_width'} * $i / $frame_opt{'scaleBar_dvdN'} ; 
		my $pY = $frame_opt{'scaleBar_y0'}; 
		my $pCol = &plotSunhh::cnvt_to_rgb( $blk_opt{'lks_col_min'} , $blk_opt{'lks_col_max'} , ( $blk_opt{'lks_col_max'}-$blk_opt{'lks_col_min'} ) * $i / $frame_opt{'scaleBar_dvdN'} , $blk_opt{'lks_col_rgb'} ); 
		$svg->rectangle(
			'x' => $pX, 'y' => $pY, 
			'width'  => $stepX, 
			'height' => $frame_opt{'scaleBar_height'}, 
			'stroke' => 'none', 
			'stroke-width' => 0.1 , 
			'fill'   => $pCol, 
		); 
	}
	$svg->text(
		'x' => $frame_opt{'scaleBar_x0'}-1, 
		'y' => $frame_opt{'scaleBar_y0'}, 
		-cdata      => "$blk_opt{'lks_col_min'}", 
		'font-size' => $frame_opt{'scaleBar_height'} * 0.8, 
		'stroke'    => 'none', 
		'fill'      => join( '', "rgb(", join(',', @{$blk_opt{'lks_col_rgb'}[0]}),")" ) ,
		'alignment-baseline' => 'text-before-edge', 
		'text-anchor'        => 'end', 
	); 
	$svg->text(
		'x' => $frame_opt{'scaleBar_x0'} + $frame_opt{'scaleBar_width'}+1, 
		'y' => $frame_opt{'scaleBar_y0'}, 
		-cdata      => "$blk_opt{'lks_col_max'}", 
		'font-size' => $frame_opt{'scaleBar_height'} * 0.8, 
		'stroke'    => 'none', 
		'fill'      => join( '', "rgb(", join(',', @{$blk_opt{'lks_col_rgb'}[-1]}),")" ) ,
		'alignment-baseline' => 'text-before-edge', 
		'text-anchor'        => 'start', 
	); 
}


##############################
# s02 - Plot small boxex for each chromosome pair. Later for label drawing. 
##############################
for ( my $xi = 0; $xi < @{$chrLisX{'arr'}} ; $xi ++ ) {
	my ($chrX_ID, $chrX_len, $chrX_cumL, $chrX_rI) = @{$chrLisX{'arr'}[$xi]}; 
	my $sP_x = (&mcsSunhh::chrP_to_plotP(
		'chrLis' => \%chrLisX , 'chrID' => $chrX_ID , 'chrP' => 0, 
		'beginPlotP' => 0, 
		'BpPoint'    => $frame_opt{'xAxis_BpPoint'} , 
		'SepBorder'  => $frame_opt{'xAxis_SepBorder'} , 
	))[ $chrX_rI ]; 
	$sP_x = $frame_lx + $sP_x ; 
	my $eP_x = (&mcsSunhh::chrP_to_plotP(
		'chrLis' => \%chrLisX , 'chrID' => $chrX_ID , 'chrP' => $chrX_len, 
		'beginPlotP' => 0, 
		'BpPoint'    => $frame_opt{'xAxis_BpPoint'} , 
		'SepBorder'  => $frame_opt{'xAxis_SepBorder'} , 
	))[ $chrX_rI ]; 
	$eP_x = $frame_lx + $eP_x ; 
	for ( my $yi = 0; $yi < @{$chrLisY{'arr'}} ; $yi ++ ) {
		my ($chrY_ID, $chrY_len, $chrY_cumL, $chrY_rI) = @{$chrLisY{'arr'}[$yi]}; 
		my $sP_y = (&mcsSunhh::chrP_to_plotP(
			'chrLis' => \%chrLisY , 'chrID' => $chrY_ID , 'chrP' => 0, 
			'beginPlotP' => 0, 
			'BpPoint'    => $frame_opt{'yAxis_BpPoint'} , 
			'SepBorder'  => $frame_opt{'yAxis_SepBorder'} , 
		))[ $chrY_rI ]; 
		$sP_y = $frame_by - $sP_y ; 
		my $eP_y = (&mcsSunhh::chrP_to_plotP(
			'chrLis' => \%chrLisY , 'chrID' => $chrY_ID , 'chrP' => $chrY_len, 
			'beginPlotP' => 0, 
			'BpPoint'    => $frame_opt{'yAxis_BpPoint'} , 
			'SepBorder'  => $frame_opt{'yAxis_SepBorder'} , 
		))[ $chrY_rI ]; 
		$eP_y = $frame_by - $eP_y ; 
		#### s02.01 - plot borders
		if ( defined $grps{'xAxis_Border'} ) {
			$grps{'xAxis_Border'}->line(
				'id' => "xAxis_Border_s:$xi:$yi", 
				'x1' => $sP_x , 'y1' => $sP_y , 
				'x2' => $sP_x , 'y2' => $eP_y , 
			); 
			$grps{'xAxis_Border'}->line(
				'id' => "xAxis_Border_e:$xi:$yi", 
				'x1' => $eP_x , 'y1' => $sP_y , 
				'x2' => $eP_x , 'y2' => $eP_y , 
			); 
		}
		if ( defined $grps{'yAxis_Border'} ) {
			$grps{'yAxis_Border'}->line(
				'id' => "yAxis_Border_s:$xi:$yi", 
				'x1' => $sP_x , 'y1' => $sP_y , 
				'x2' => $eP_x , 'y2' => $sP_y , 
			); 
			$grps{'yAxis_Border'}->line(
				'id' => "yAxis_Border_e:$xi:$yi", 
				'x1' => $sP_x , 'y1' => $eP_y , 
				'x2' => $eP_x , 'y2' => $eP_y , 
			); 
		}
		#### s02.02 - plot ticks : 
		if ( $yi == 0 and defined $grps{'xAxis_Tick'} ) {
			for ( my $i=0; $i<=$chrX_len; $i+=$frame_opt{'xAxis_TickStep'} ) {
				my $show_v = $i/$frame_opt{'xAxis_TickUnitLen'}; 
				$show_v .= "$frame_opt{'xAxis_TickUnit'}"; 
				$grps{'xAxis_Tick'}->line(
					'id' => "xAxis_tickBar:$chrX_ID:$chrX_rI:$i", 
					'x1' => $i / $frame_opt{'xAxis_BpPoint'}, 'y1' => 0, 
					'x2' => $i / $frame_opt{'xAxis_BpPoint'}, 'y2' => $frame_opt{'xAxis_TickLen'}, 
					'transform' => join('', "translate($sP_x $sP_y)"), 
				); 
				$grps{'xAxis_Tick'}->text(
					'id' => "xAxis_tickTxt:$chrX_ID:$chrX_rI:$i", 
					'x'  => $i / $frame_opt{'xAxis_BpPoint'}, 'y'  => $frame_opt{'xAxis_TickLen'}, 
					-cdata => $show_v, 
					'transform' => join('', "translate($sP_x $sP_y)"), 
					'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
					'stroke' => 'none', 
				); # 'alignment-baseline' doesn't work in groups. 
			}
		}
		if ( $xi == 0 and defined $grps{'yAxis_Tick'} ) {
			for ( my $i=0; $i<=$chrY_len; $i+=$frame_opt{'yAxis_TickStep'} ) {
				my $show_v = $i/$frame_opt{'yAxis_TickUnitLen'}; 
				$show_v .= "$frame_opt{'yAxis_TickUnit'}"; 
				$grps{'yAxis_Tick'}->line(
					'id' => "yAxis_tickBar:$chrY_ID:$chrY_rI:$i", 
					'x1' => 0,                            'y1' => -1 * $i / $frame_opt{'yAxis_BpPoint'}, 
					'x2' => -$frame_opt{'yAxis_TickLen'}, 'y2' => -1 * $i / $frame_opt{'yAxis_BpPoint'}, 
					'transform' => join('', "translate($sP_x $sP_y)"), 
				); 
				my $txt_xV = -$frame_opt{'yAxis_TickLen'}; 
				my $txt_yV = -1 * $i / $frame_opt{'yAxis_BpPoint'}; 
				$grps{'yAxis_Tick'}->text(
					'id' => "yAxis_tickTxt:$chrY_ID:$chrY_rI:$i", 
					'x'  => $txt_xV , 'y' => $txt_yV , 
					-cdata => $show_v, 
					'transform' => join('', "translate($sP_x $sP_y) rotate(90 $txt_xV $txt_yV)"), 
					'alignment-baseline'=> $frame_opt{'yAxis_TxtVertAln'}, 
					'stroke' => 'none', 
				); # 'alignment-baseline' doesn't work in groups. 
			}
		}
		#### s02.03 - Plot chromosome ID text 
		if ( $yi == 0 and $frame_opt{'xAxis_ChrIDFontSize'} > 0 ) {
			$grps{'xAxis_Label'}->text(
				'x'    => ($sP_x + $eP_x)/2 , 
				'y'    => $frame_by + $frame_opt{'xAxis_TickLen'} + $frame_opt{'xAxis_TickFontSize'} + $frame_opt{'xAxis_ChrIDFontSize'} * 0.3 , 
				-cdata => "$chrX_ID" , 
				'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
				'font-size' => $frame_opt{'xAxis_ChrIDFontSize'} , 
			); 
		}
		if ( $xi == 0 and $frame_opt{'yAxis_ChrIDFontSize'} > 0 ) {
			my $chrID_xV = $frame_lx - $frame_opt{'yAxis_TickLen'} - $frame_opt{'yAxis_TickFontSize'} - $frame_opt{'yAxis_ChrIDFontSize'} * 0.3; 
			my $chrID_yV = ($sP_y+$eP_y)/2 ; 
			$grps{'yAxis_Label'}->text(
				'x' => $chrID_xV, 'y' => $chrID_yV, 
				-cdata=>"$chrY_ID", 
				'transform' => join('', "rotate(90 $chrID_xV $chrID_yV)"), 
				'alignment-baseline'=> $frame_opt{'yAxis_TxtVertAln'}, 
				'font-size' => $frame_opt{'yAxis_ChrIDFontSize'} , 
			); # For chromosome ID text 
		}
		#### s02.04 - Plot syntenic blocks
		# $opts{'show_time'} and &tsmsg("[Msg] Plotting blocks for pair : $chrX_ID ($xi) - $chrY_ID ($yi)\n"); 
		for my $blkR (@tabInfo) {
			my ($blkID, $chrLoc_1, $chrLoc_2, $strand, $alnScore, $alnEvalue, $alnNumber, $geneLis_1, $geneLis_2, $kaLis, $ksLis, $kaksLis, @rest) = @$blkR; 
			my ($loc1_ID, $loc1_S, $loc1_E) = @$chrLoc_1; 
			my ($loc2_ID, $loc2_S, $loc2_E) = @$chrLoc_2; 
			$strand eq '-' and ($loc2_S, $loc2_E) = ($loc2_E, $loc2_S); 
			if ( $chrX_ID eq $loc1_ID and $chrY_ID eq $loc2_ID ) { 
				my ( $rel_sPx, $rel_ePx ) = ( $loc1_S / $frame_opt{'xAxis_BpPoint'}, $loc1_E / $frame_opt{'xAxis_BpPoint'} ); 
				my ( $rel_sPy, $rel_ePy ) = ( -1 * $loc2_S / $frame_opt{'yAxis_BpPoint'}, -1 * $loc2_E / $frame_opt{'yAxis_BpPoint'} ); 
				$rest[0] =~ m/^(NA|NAN)$/i and $rest[0] = $blk_opt{'lks_col_max'}; 
				my $link_col = &plotSunhh::cnvt_to_rgb( $blk_opt{'lks_col_min'} , $blk_opt{'lks_col_max'} , $rest[0] , $blk_opt{'lks_col_rgb'} ); 
				$grps{'blk_line'}->line(
					'id' => "blk_line:blkID_${blkID}:$xi:$yi:1to2", 
					'x1' => $rel_sPx , 'y1' => $rel_sPy , 
					'x2' => $rel_ePx , 'y2' => $rel_ePy , 
					'transform' => join('', "translate($sP_x, $sP_y)"), 
					'stroke'    => $link_col, 
				); 
				for ( my $j=0; $j<@$geneLis_1; $j++ ) {
					my @gx_loc = @{ $gen2loc{$geneLis_1->[$j]} }; 
					my @gy_loc = @{ $gen2loc{$geneLis_2->[$j]} }; 
					my ( $rel_pPx, $rel_pPy ) = ( ($gx_loc[1]+$gx_loc[2])/2 / $frame_opt{'xAxis_BpPoint'} , -1 * ($gy_loc[1]+$gy_loc[2])/2 / $frame_opt{'yAxis_BpPoint'} ); 
					$ksLis->[$j] =~ m/^(NA|NAN)$/i and $ksLis->[$j] = $blk_opt{'pks_col_max'}; 
					my $point_col = &plotSunhh::cnvt_to_rgb( $blk_opt{'pks_col_min'} , $blk_opt{'pks_col_max'} , $ksLis->[$j], $blk_opt{'pks_col_rgb'} ); 
					$grps{'blk_point'}->circle(
						'id' => "blk_point:blkID_${blkID}:$xi:$yi:1to2:Gene_${j}" , 
						'cx' => $rel_pPx, 'cy' => $rel_pPy, 
						'r'  => $blk_opt{'p_radius'}, 
						'transform' => join('', "translate($sP_x, $sP_y)") , 
						'stroke' => $point_col, 
						'fill'   => $point_col, 
					); 
				}
			}# if ( chrX_ID eq loc1_ID 
			if ( $chrX_ID eq $loc2_ID and $chrY_ID eq $loc1_ID ) {
				my ( $rel_sPx, $rel_ePx ) = ( $loc2_S / $frame_opt{'xAxis_BpPoint'}, $loc2_E / $frame_opt{'xAxis_BpPoint'} ); 
				my ( $rel_sPy, $rel_ePy ) = ( -1 * $loc1_S / $frame_opt{'yAxis_BpPoint'}, -1 * $loc1_E / $frame_opt{'yAxis_BpPoint'} ); 
				$rest[0] =~ m/^(NA|NAN)$/i and $rest[0] = $blk_opt{'lks_col_max'}; 
				my $link_col = &plotSunhh::cnvt_to_rgb( $blk_opt{'lks_col_min'} , $blk_opt{'lks_col_max'} , $rest[0] , $blk_opt{'lks_col_rgb'} ); 
				$grps{'blk_line'}->line(
					'id' => "blk_line:blkID_${blkID}:$xi:$yi:2to1", 
					'x1' => $rel_sPx , 'y1' => $rel_sPy , 
					'x2' => $rel_ePx , 'y2' => $rel_ePy , 
					'transform' => join('', "translate($sP_x, $sP_y)"), 
					'stroke'    => $link_col, 
				); 
				for ( my $j=0; $j<@$geneLis_1; $j++ ) {
					my @gx_loc = @{ $gen2loc{$geneLis_2->[$j]} }; 
					my @gy_loc = @{ $gen2loc{$geneLis_1->[$j]} }; 
					my ( $rel_pPx, $rel_pPy ) = ( ($gx_loc[1]+$gx_loc[2])/2 / $frame_opt{'xAxis_BpPoint'} , -1 * ($gy_loc[1]+$gy_loc[2])/2 / $frame_opt{'yAxis_BpPoint'} ); 
					$ksLis->[$j] =~ m/^(NA|NAN)$/i and $ksLis->[$j] = $blk_opt{'pks_col_max'}; 
					my $point_col = &plotSunhh::cnvt_to_rgb( $blk_opt{'pks_col_min'} , $blk_opt{'pks_col_max'} , $ksLis->[$j], $blk_opt{'pks_col_rgb'} ); 
					$grps{'blk_point'}->circle(
						'id' => "blk_point:blkID_${blkID}:$xi:$yi:2to1:Gene_${j}" , 
						'cx' => $rel_pPx, 'cy' => $rel_pPy, 
						'r'  => $blk_opt{'p_radius'}, 
						'transform' => join('', "translate($sP_x, $sP_y)") , 
						'stroke' => $point_col, 
						'fill'   => $point_col, 
					); 
				}
			}# if ( chrX_ID eq loc2_ID 
		}
		#### s02.05 - Plot shadow lines. 
		if ( $yi == 0 and @{$shad_opt{'xAxis_shadow_setopt'}} > 0 ) {
			my $hori_yV; 
			for ( my $j=0; $j<@{$shad_opt{'xAxis_shadow_setopt'}}; $j++ ) {
				my ( $shad_ID, $shad_fn, $shad_col, $shad_height, $shad_locHR ) = @{ $shad_opt{'xAxis_shadow_setopt'}[$j] }; 
				my $shad_yV = $frame_by + $frame_opt{'xAxis_TickLen'} + $frame_opt{'xAxis_TickFontSize'} + $frame_opt{'xAxis_ChrIDFontSize'}*2 + $frame_opt{'xAxis_LabFontSize'}*1.5 ; 
				$hori_yV //= $shad_yV ; 
				$shad_yV = $hori_yV; 
				$svg->rectangle(
					'x' => $sP_x , 'y' => $shad_yV, 
					'width'  => $eP_x - $sP_x , 
					'height' => $shad_height , 
					'stroke'       => $other_opt{'shadow_stroke'} , 
					'stroke-width' => $other_opt{'shadow_strokeWid'}, 
					'fill'         => $other_opt{'shadow_fill'}, 
				); 
				if ( $xi == $#{$chrLisX{'arr'}} ) {
					$svg->text(
						'x' => $eP_x + 3 + $shad_height * 0.8 , 'y' => $shad_yV , 
						-cdata => "$shad_ID" , 
						'fill'               => $shad_col, 
						'stroke'             => 'none' , 
						'font-size'          => $shad_height * 0.8, 
						'text-anchor'        => 'start', 
						'alignment-baseline' => 'text-before-edge', 
					); 
				}
				defined $shad_locHR->{$chrX_ID} or next; 
				for my $t_se ( @{$shad_locHR->{$chrX_ID}} ) {
					my $rel_sP = $t_se->[0] / $frame_opt{'xAxis_BpPoint'} ; 
					my $rel_eP = $t_se->[1] / $frame_opt{'xAxis_BpPoint'} ; 
					my $shad_xV = $sP_x + $rel_sP; 
					$svg->rectangle(
						'x' => $shad_xV, 'y' => $shad_yV , 
						'width'  => $rel_eP - $rel_sP, 
						'height' => $shad_height, 
						'stroke'       => 'none', 
						'stroke-width' => 0, 
						'fill'         => $shad_col, 
					); 
				}
				$hori_yV += ( $shad_height * 1.2 ) ; 
			}
		}
		if ( $xi == 0 and @{$shad_opt{'yAxis_shadow_setopt'}} > 0 ) {
			my $vert_xV; 
			for ( my $j=0; $j<@{$shad_opt{'yAxis_shadow_setopt'}}; $j++ ) {
				my ( $shad_ID, $shad_fn, $shad_col, $shad_height, $shad_locHR ) = @{ $shad_opt{'yAxis_shadow_setopt'}[$j] }; 
				my $shad_xV = $frame_lx - $frame_opt{'yAxis_TickLen'} - $frame_opt{'yAxis_TickFontSize'} - $frame_opt{'yAxis_ChrIDFontSize'}*2 - $frame_opt{'yAxis_LabFontSize'}*1.5 ; 
				$vert_xV //= $shad_xV ; 
				$shad_xV = $vert_xV; 
				$svg->rectangle(
					'x' => $shad_xV - $shad_height , 'y' => $eP_y, 
					'width'  => $shad_height , 
					'height' => $sP_y - $eP_y , 
					'stroke'       => $other_opt{'shadow_stroke'} , 
					'stroke-width' => $other_opt{'shadow_strokeWid'}, 
					'fill'         => $other_opt{'shadow_fill'}, 
				); 
				if ( $yi == 0 ) {
					my $tmp_yV = $sP_y + 3 + $shad_height * 0.8 ; 
					$svg->text(
						'x' => $shad_xV , 'y' => $tmp_yV , 
						-cdata => "$shad_ID" , 
						'fill'               => $shad_col, 
						'stroke'             => 'none' , 
						'font-size'          => $shad_height * 0.8, 
						'text-anchor'        => 'start', 
						'alignment-baseline' => 'text-before-edge', 
						'transform'          => "rotate(90 $shad_xV $tmp_yV)", 
					); 
				}
				defined $shad_locHR->{$chrY_ID} or next; 
				for my $t_se ( @{$shad_locHR->{$chrY_ID}} ) {
					my $rel_sP = $t_se->[0] / $frame_opt{'yAxis_BpPoint'} ; 
					my $rel_eP = $t_se->[1] / $frame_opt{'yAxis_BpPoint'} ; 
					my $shad_yV = $sP_y - $rel_eP; 
					$svg->rectangle(
						'x' => $shad_xV - $shad_height, 'y' => $shad_yV , 
						'width'  => $shad_height , 
						'height' => $rel_eP - $rel_sP , 
						'stroke'       => 'none', 
						'stroke-width' => 0, 
						'fill'         => $shad_col, 
					); 
				}
				$vert_xV -= ( $shad_height * 1.2 ) ; 
			}
		}
		#### s02.06 - Plot additional lines. (From -add_lineX or -add_lineY)
		if ( $yi == 0 and defined $addLineX{$chrX_ID} ) {
			for my $cp ( @{$addLineX{$chrX_ID}} ) {
				my $rel_P = $cp / $frame_opt{'xAxis_BpPoint'}; 
				$grps{'xAxis_Border'}->line(
					'x1' => $sP_x + $rel_P , 'y1' => $frame_by , 
					'x2' => $sP_x + $rel_P , 'y2' => $frame_ty , 
				); 
			}
		}
		if ( $xi == 0 and defined $addLineY{$chrY_ID} ) {
			for my $cp ( @{$addLineY{$chrY_ID}} ) {
				my $rel_P = $cp / $frame_opt{'yAxis_BpPoint'}; 
				$grps{'yAxis_Border'}->line(
					'x1' => $frame_lx , 'y1' => $sP_y - $rel_P , 
					'x2' => $frame_rx , 'y2' => $sP_y - $rel_P , 
				); 
			}
		}
	}# for my $chrY_ar
}# for my $chrX_ar

##############################
# s03 - Plot X-Y axes labels. 
##############################
if ( $img_opt{'xAxisLabel'} ne '' ) {
	$grps{'xAxis_Label'}->text(
		'x' => ($frame_lx + $frame_rx) / 2, 'y' => $frame_by + $frame_opt{'xAxis_TickLen'} + $frame_opt{'xAxis_TickFontSize'} + $frame_opt{'xAxis_ChrIDFontSize'}*2 , 
		-cdata => $img_opt{'xAxisLabel'}, 
		'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
	); 
}# Label
if ( $img_opt{'yAxisLabel'} ne '' ) {
	my $xV = $frame_lx - $frame_opt{'yAxis_TickLen'} - $frame_opt{'yAxis_TickFontSize'} - $frame_opt{'yAxis_ChrIDFontSize'}*2 ; 
	my $yV = ($frame_by + $frame_ty) / 2; 
	$grps{'yAxis_Label'}->text(
		'x' => $xV, 'y' => $yV, 
		-cdata => $img_opt{'yAxisLabel'}, 
		'transform' => join('', "rotate(90 $xV $yV)"), 
		'alignment-baseline'=> $frame_opt{'yAxis_TxtVertAln'}, 
	); 
}# Label




##############################
# sXX - Print out 
##############################
print {$outSvgFh} $svg->xmlify; 

$opts{'show_time'} and &tsmsg("[Msg] Script [$0] all done.\n"); 

################################################################################
#    Sub-routines 
################################################################################

################################################################################
#    Inner sub-routines 
################################################################################
=head1 _pp2xy( 'Px'=>[chrID_x, pos_x, \%chrLis_x] , 'Py'=>[chrID_y, pos_y, \%chrLis_y] )

Return       : ( [xV_1, yV_1], [xV_2, yV_2], ... )

=cut
sub _pp2xy {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	$parm{'beginPlotP'} //= 0; 
	my %xPar = %parm; 
	my %yPar = %parm; 
	for my $t (qw/Px Py/) { delete($xPar{$t}); delete($yPar{$t}); }

	$xPar{'chrID'}  = $parm{'Px'}[0]; 
	$xPar{'chrP'}   = $parm{'Px'}[1]; 
	$xPar{'chrLis'} = $parm{'Px'}[2]; 
	$yPar{'chrID'}  = $parm{'Py'}[0]; 
	$yPar{'chrP'}   = $parm{'Py'}[1]; 
	$yPar{'chrLis'} = $parm{'Py'}[2]; 
	$xPar{'beginPlotP'} //= 0; 
	$xPar{'BpPoint'} //= $frame_opt{'xAxis_BpPoint'}; 
	$yPar{'beginPlotP'} //= 0; 
	$xPar{'SepBorder'} //= $frame_opt{'xAxis_SepBorder'}; 
	$yPar{'BpPoint'} //= $frame_opt{'yAxis_BpPoint'}; 
	$yPar{'SepBorder'} //= $frame_opt{'yAxis_SepBorder'}; 
	
	my @xP = &mcsSunhh::chrP_to_plotP(%xPar); 
	my @yP = &mcsSunhh::chrP_to_plotP(%yPar); 
	my @back; 
	for my $tx (@xP) {
		for my $ty (@yP) {
			push(@back, [$tx, $ty]); 
		}
	}
	
	return(@back); 
}# _pp2xy 

=head1 _readInAddLine ( $file_name ) 

Input        : Format : chrID \\t position \\n

Return       : ( { chrID => [pos_1, pos_2], ... } )

=cut
sub _readInAddLine {
	my $fn = shift; 
	my %back; 
	defined $fn or return (\%back); 
	my $fh = &openFH( $fn, '<' ); 
	while (&wantLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		push(@{$back{$ta[0]}}, $ta[1]); 
	}
	close($fh); 
	return(\%back); 
} # _readInAddLine() 

sub _shadow_setopt {
	for my $tk (qw/xAxis_shadow_setopt yAxis_shadow_setopt/) {
		defined $opts{$tk} or do { $shad_opt{$tk} = []; next; }; 
		for my $t1 (@{$opts{$tk}}) {
			for my $t2 ( split(/;/, $t1) ) {
				$t2 =~ m/^(\S+)=(.*)$/i or &stopErr("[Err] Failed at [$t2]\n"); 
				my ($k, $v) = ($1, $2); 
				(defined $v and $v ne '') or next; 
				my @ta = split(/:/, $v); 
				(defined $ta[1] and $ta[1] ne '') or $ta[1] = 'black'; 
				(defined $ta[2] and $ta[2] ne '') or $ta[2] = 10; 
				my %blk; 
				my $fh = &openFH( $ta[0], '<' ); 
				# shadow_01 file format : ChromeID \\t ChromeStart \\t ChromeEnd \\n 
				while (&wantLineC($fh)) {
					my @tb = &splitL("\t", $_); 
					$tb[1] > $tb[2] and @tb[1,2] = @tb[2,1]; 
					push(@{$blk{$tb[0]}}, [@tb[1,2]]); 
				}
				close($fh); 
				push(@{$shad_opt{$tk}}, [$k, @ta, \%blk]); 
			}
		}
	}
	return; 
} # _shadow_setopt() 

sub _frame_setopt {
	&_setoptAR( \%frame_opt, $opts{'frame_setopt'}, 1 ); 
	&_setoptAR( \%frame_opt, $default_string{'frame_setopt'}, 0); 
	unless ( defined $frame_opt{'xAxis_TickUnitLen'} ) {
		if ( $frame_opt{'xAxis_TickUnit'} eq '' ) {
			$frame_opt{'xAxis_TickUnitLen'} = 1; 
		} elsif ( $frame_opt{'xAxis_TickUnit'} =~ m/^k$/i ) {
			$frame_opt{'xAxis_TickUnitLen'} = 1e3; 
		} elsif ( $frame_opt{'xAxis_TickUnit'} =~ m/^M$/i ) {
			$frame_opt{'xAxis_TickUnitLen'} = 1e6; 
		} else {
			&stopErr("[Err] Unknown xAxis_TickUnit [$frame_opt{'xAxis_TickUnit'}]\n"); 
		}
	}
	unless ( defined $frame_opt{'yAxis_TickUnitLen'} ) {
		if ( $frame_opt{'yAxis_TickUnit'} eq '' ) {
			$frame_opt{'yAxis_TickUnitLen'} = 1; 
		} elsif ( $frame_opt{'yAxis_TickUnit'} =~ m/^k$/i ) {
			$frame_opt{'yAxis_TickUnitLen'} = 1e3; 
		} elsif ( $frame_opt{'yAxis_TickUnit'} =~ m/^M$/i ) {
			$frame_opt{'yAxis_TickUnitLen'} = 1e6; 
		} else {
			&stopErr("[Err] Unknown yAxis_TickUnit [$frame_opt{'yAxis_TickUnit'}]\n"); 
		}
	}
	return; 
}# _frame_setopt() 

sub _blk_setopt {
	&_setoptAR( \%blk_opt, $opts{'blk_setopt'}, 1 ); 
	&_setoptAR( \%blk_opt, $default_string{'blk_setopt'}, 0); 

	for my $t (qw/l p/) {
		my $c1 = "${t}ks_col_rgb"; 
		my $c2 = "${t}_col"; 
		( defined $blk_opt{$c1} and $blk_opt{$c1} ne '' ) or $blk_opt{$c1} = $blk_opt{$c2}; 
		my $t1 = $blk_opt{$c1}; 
		$blk_opt{$c1} = []; 
		for my $t2 (split(/:/, $t1)) {
			my $rgb ; 
			if ( $t2 =~ m/^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i ) {
				push(@{$blk_opt{$c1}}, [$1,$2,$3]); 
			} elsif ( $rgb = &plotSunhh::_rgb_color( 'in_color' => "html_name:$t2" , 'out_fmt' => 'rgb' ) ) {
				$rgb =~ m/^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i or &stopErr("[Err] Bad RGB format [$rgb]\n"); 
				push(@{$blk_opt{$c1}}, [$1,$2,$3]); 
			} elsif ( $rgb = &plotSunhh::_rgb_color( 'in_color' => "hex_name:$t2" , 'out_fmt' => 'rgb' ) ) {
				$rgb =~ m/^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i or &stopErr("[Err] Bad RGB format [$rgb]\n"); 
				push(@{$blk_opt{$c1}}, [$1,$2,$3]); 
			} else {
				&stopErr("[Err] Failed to parse color [$t2]\n"); 
			}
		}
	}
	
	return; 
}# _blk_setopt() 


sub _img_setopt {
	&_setoptAR( \%img_opt, $opts{'img_setopt'}, 1 ); 
	&_setoptAR( \%img_opt, $default_string{'img_setopt'}, 0 ); 
	return; 
}# _img_setopt() 

sub _other_setopt {
	&_setoptAR( \%other_opt, $opts{'other_setopt'}, 1 ); 
	&_setoptAR( \%other_opt, $default_string{'other_setopt'}, 0 ); 
}# _other_setopt() 

sub _setoptAR {
	my ($hr, $txtAR, $overwrite) = @_; 
	$overwrite //= 1; 
	defined $txtAR or return; 
	for my $t1 (@$txtAR) {
		defined $t1 or next; 
		&_setopt($hr, $t1, $overwrite); 
	}
	return; 
} # _setoptAR()

sub _setopt {
	my ($hr, $txt, $overwrite) = @_; 
	$overwrite //= 1; 
	defined $txt or return; 
	for my $t2 (split(";", $txt)) {
		$t2 =~ s!^\s+|\s+$!!g; 
		$t2 eq '' and next; 
		$t2 =~ m/^(\S+)=(.*)$/i or &stopErr("[Err] Failed to parse opt [$t2]\n"); 
		my ($k, $v) = ($1, $2); 
		$v //= ''; 
		( $overwrite or !(defined $hr->{$k}) ) and $hr->{$k} = $v; 
	}
	return; 
}# _setopt



=head1 _seg2xy ( 'Ps' => [chrID_x_1, pos_x_1, chrID_y_1, pos_y_1] , 'Pe' => [chrID_x_2, pos_x_2, chrID_y_2, pos_y_2], 'chrLisX' => \%chrLisX, 'chrLisX' => \%chrLisY ) 

Return       : ( [sVx_1, sVy_1, eVx_1, eVy_1], [sVx_2, sVy_2, eVx_2, eVy_2], ... )

=cut
sub _seg2xy {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	my @back; 
	my %parm_Ps = %parm; 
	$parm_Ps{'Px'} = [ $parm_Ps{'Ps'}[0], $parm_Ps{'Ps'}[1], $parm_Ps{'chrLisX'} ]; 
	$parm_Ps{'Py'} = [ $parm_Ps{'Ps'}[2], $parm_Ps{'Ps'}[3], $parm_Ps{'chrLisY'} ]; 
	my @xy_Ps = &_pp2xy( %parm_Ps ); 
	my %parm_Pe = %parm; 
	$parm_Pe{'Px'} = [ $parm_Pe{'Pe'}[0], $parm_Pe{'Pe'}[1], $parm_Pe{'chrLisX'} ]; 
	$parm_Pe{'Py'} = [ $parm_Pe{'Pe'}[2], $parm_Pe{'Pe'}[3], $parm_Pe{'chrLisY'} ]; 
	my @xy_Pe = &_pp2xy( %parm_Pe ); 
	for my $ts (@xy_Ps) {
		for my $te (@xy_Pe) {
			push(@back, [ $ts->[0], $ts->[1], $te->[0], $te->[1] ]); 
		}
	}
	
	return(@back); 
}# _seg2xy() 

