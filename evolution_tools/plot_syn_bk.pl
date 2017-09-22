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
	"out_svg:s", # 
	# Input alignments : 
	"in_tab:s", # result from "follow_mcscan.pl -aln2table". 
	# Input gene positions : MCScan formatted .gff . 
	"in_pos:s", # Format : Chr_ID  \\t  Gene_ID  \\t  Gene_Start  \\t  Gene_End (unstranded)
	# Input X-Y axis
	"chrLis_x:s", "chrLis_y:s", # Format : Chr_ID \\t Chr_Len \\n
	# For frame options. 
	"frame_setopt:s@", # 
	# For syntenic blocks 
	"blk_setopt:s@", # Can be assigned multiple times. 
	# For shadow lines : 
	"xAxis_shadow_setopt:s@", "yAxis_shadow_setopt:s@", 
	# For additional lines
	"add_lineX:s", "add_lineY:s", # Add border lines for x/y axises. format : Chr_ID \\t Chr_Position \\n ; The border line will be drawn at Chr_ID - Chr_Position ; 
	# Other settings 
	"img_setopt:s@", # Can be assigned multiple times. $default_string{'img_setopt'}
); 

################################################################################
#    Setup basic parameters. 
################################################################################
my $outSvgFh = \*STDOUT ; 
defined $opts{'out_svg'} and $outSvgFh = &openFH( $opts{'out_svg'}, '>' ); 


my %default_string; 
$default_string{'img_setopt'} = ['title=;horizMargin=50;vertMargin=50;xAxisLabel=Genome 1;yAxisLabel=Genome 2']; 
$default_string{'blk_setopt'} = [
	'l_lwd=1;l_col=black;',                                                # lines linking points within block 
	'lks_col_min=0;lks_col_max=3;lks_col_rgb=blue:cyan:green:yellow:red;', # block's average Ks : ''
	'p_radius=1;p_lwd=0;',                                                 # points standing for gene pairs : ''
	'pks_col_min=0;pks_col_max=3;pks_col_rgb=blue:cyan:green:yellow:red;', # gene pairs' Ks     : ''
]; 
$default_string{'frame_setopt'} = [
	'title_FontSize=40;title_FontFam=ArialNarrow;title_FontWeight=bold;title_HoriAln=middle;title_VerAln=text-after-edge;', 
	
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

$opts{'help'} and &usage(); 
defined $opts{'in_tab'} or &usage(); 

sub usage {
	print <<HH; 
################################################################################
#  perl $0 
#   -help
#
#  Required : 
#   -in_tab      out.collinearity.tab 
#   -in_pos      mcscanX.gff . Format : Chr_ID  \\t  Gene_ID  \\t  Gene_Start  \\t  Gene_End (unstranded)
#   -chrLis_x    Format : Chr_ID \\t Chr_Len \\n
#   -chrLis_y    Format : Chr_ID \\t Chr_Len \\n
#
#  .....
################################################################################
HH
	exit(1); 
}


################################################################################
#    Read in data files. 
################################################################################
my @tabInfo = @{ &mcsSunhh::_readInAlnTbl( $opts{'in_tab'} ) }; 
my @tabHeader; $tabInfo[0] =~ m/^BlkID$/i and do { @tabHeader=@{$tabInfo[0]}; shift(@tabInfo); }; 
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
my $svg_width = &mathSunhh::max( 
    &mcsSunhh::chrP_to_plotP( 
      'chrLis' => \%chrLisX, 'chrID'=>$chrLisX{'arr'}[-1][0], 'chrP'=>$chrLisX{'arr'}[-1][1], 
      'beginPlotP' => $img_opt{'horizMargin'}, 
      'BpPoint' => $frame_opt{'xAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'xAxis_SepBorder'}, 
    ) 
  ) 
  + $img_opt{'horizMargin'}
; 
my $svg_height = &mathSunhh::max( 
    &mcsSunhh::chrP_to_plotP( 
      'chrLis' => \%chrLisY, 'chrID'=>$chrLisY{'arr'}[-1][0], 'chrP'=>$chrLisY{'arr'}[-1][1], 
      'beginPlotP' => $img_opt{'vertMargin'}, 
      'BpPoint' => $frame_opt{'yAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'yAxis_SepBorder'}, 
    )
  ) 
  + $img_opt{'vertMargin'}
; 
my $svg = SVG->new( 'width' => $svg_width, 'height' => $svg_height ); 

##############################
# SVG groups : 
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
); 
$grps{'blk_point'} = $svg->group(
	'id'                => "blk_point", 
	'stroke-width'      => $blk_opt{'p_lwd'}, 
	'stroke'            => 'black', 
	'opacity'           => 1, 
); 

##############################
# s01 - Main frame
##############################
my $frame_lx = $img_opt{'horizMargin'}; # left-X 
my $frame_by = $svg_height - $img_opt{'vertMargin'}; # bottom-Y
my $frame_rx = $svg_width - $img_opt{'horizMargin'}; # right-X
my $frame_ty = $img_opt{'vertMargin'}; # top-Y
my $frame_width = 
$svg->rectangle(
	'x'      => $frame_lx, 'y' => $frame_ty, 
	'width'  => $frame_rx - $frame_lx, 
	'height' => $frame_by - $frame_ty, 
	'id'     => "Main_Frame", 
	'stroke' => 'black', 
	'fill'   => 'none', 
); 
##############################
# s02 - x axis : Label , chrIDs , borders, and ticks. 
##############################
for my $ele1 ( @{$chrLisX{'arr'}} ) {
	my ($id, $len, $cumP, $rI) = @$ele1; 
	my @sP = &mcsSunhh::chrP_to_plotP(
		'chrLis' => \%chrLisX, 'chrID'=>$id, 'chrP'=>1, 
		'beginPlotP' => 0, 
		'BpPoint' => $frame_opt{'xAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'xAxis_SepBorder'}, 
	); 
	my @eP = &mcsSunhh::chrP_to_plotP(
		'chrLis' => \%chrLisX, 'chrID'=>$id, 'chrP'=>$len, 
		'beginPlotP' => 0, 
		'BpPoint' => $frame_opt{'xAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'xAxis_SepBorder'}, 
	); 
	if ( defined $grps{'xAxis_Border'} ) {
		$grps{'xAxis_Border'}->line(
			'id' => "xAxis_border_chrS:$id:1", 
			'x1' => $frame_lx + $sP[$rI], 'y1' => $frame_ty, 
			'x2' => $frame_lx + $sP[$rI], 'y2' => $frame_by, 
		); 
		$grps{'xAxis_Border'}->line(
			'id' => "xAxis_border_chrE:$id:$len", 
			'x1' => $frame_lx + $eP[$rI], 'y1' => $frame_ty, 
			'x2' => $frame_lx + $eP[$rI], 'y2' => $frame_by, 
		); 
	}# grps_xAxis_Border : borders 
	if ( defined $grps{'xAxis_Tick'} ) {
		for ( my $i=0; $i<=$len; $i+=$frame_opt{'xAxis_TickStep'} ) {
			my $show_v = $i/$frame_opt{'xAxis_TickUnitLen'}; 
			$show_v .= "$frame_opt{'xAxis_TickUnit'}"; 
			$grps{'xAxis_Tick'}->line(
				'id' => "xAxis_tick:$rI:$id:$i", 
				'x1' => ($i-1) / $frame_opt{'xAxis_BpPoint'}, 'y1' => $frame_by, 
				'x2' => ($i-1) / $frame_opt{'xAxis_BpPoint'}, 'y2' => $frame_by+$frame_opt{'xAxis_TickLen'}, 
				'transform' => join('', "translate(" , $frame_lx + $sP[$rI] , ")"), 
			); 
			$grps{'xAxis_Tick'}->text(
				'id' => "xAxis_tickTxt:$rI:$id:$i", 
				'x' => ($i-1) / $frame_opt{'xAxis_BpPoint'}, 'y' => $frame_by+$frame_opt{'xAxis_TickLen'}, 
				-cdata => $show_v, 
				'transform' => join('', "translate(" , $frame_lx + $sP[$rI], ")"), 
				'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
				'stroke' => 'none', 
			); # 'alignment-baseline' doesn't work in groups. 
		}
	}# grps_xAxis_Tick : ticks
	$grps{'xAxis_Label'}->text(
		'x' => $frame_lx + ($sP[$rI] + $eP[$rI])/2 , 'y'=> $frame_by + $frame_opt{'xAxis_TickLen'} + $frame_opt{'xAxis_TickFontSize'} + $frame_opt{'xAxis_ChrIDFontSize'} * 0.3 , 
		-cdata=>"$id", 
		'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
		'font-size' => $frame_opt{'xAxis_ChrIDFontSize'} , 
	); # For chromosome ID text 
}
if ( $img_opt{'xAxisLabel'} ne '' ) {
	$grps{'xAxis_Label'}->text(
		'x' => ($frame_lx + $frame_rx) / 2, 'y' => $frame_by + $frame_opt{'xAxis_TickLen'} + $frame_opt{'xAxis_TickFontSize'} + $frame_opt{'xAxis_ChrIDFontSize'}*2 , 
		-cdata => $img_opt{'xAxisLabel'}, 
		'alignment-baseline'=> $frame_opt{'xAxis_TxtVertAln'}, 
	); 
}# Label

##############################
# s03 - y axis : Label , chrIDs , borders, and ticks. 
##############################
for my $ele1 ( @{$chrLisY{'arr'}} ) {
	my ($id, $len, $cumP, $rI) = @$ele1; 
	my @sP = &mcsSunhh::chrP_to_plotP(
		'chrLis' => \%chrLisY, 'chrID'=>$id, 'chrP'=>1, 
		'beginPlotP' => 0, 
		'BpPoint' => $frame_opt{'yAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'yAxis_SepBorder'}, 
	); 
	my @eP = &mcsSunhh::chrP_to_plotP(
		'chrLis' => \%chrLisY, 'chrID'=>$id, 'chrP'=>$len, 
		'beginPlotP' => 0, 
		'BpPoint' => $frame_opt{'yAxis_BpPoint'}, 'SepBorder'  => $frame_opt{'yAxis_SepBorder'}, 
	); 
	if ( defined $grps{'yAxis_Border'} ) {
		$grps{'yAxis_Border'}->line(
			'id' => "yAxis_border_chrS:$id:1", 
			'x1' => $frame_lx, 'y1' => $frame_by - $sP[$rI], 
			'x2' => $frame_rx, 'y2' => $frame_by - $sP[$rI], 
		); 
		$grps{'yAxis_Border'}->line(
			'id' => "yAxis_border_chrE:$id:$len", 
			'x1' => $frame_lx, 'y1' => $frame_by - $eP[$rI], 
			'x2' => $frame_lx, 'y2' => $frame_by - $eP[$rI], 
		); 
	}# grps_yAxis_Border : borders 
	if ( defined $grps{'yAxis_Tick'} ) {
		for ( my $i=0; $i<=$len; $i+=$frame_opt{'yAxis_TickStep'} ) {
			my $show_v = $i/$frame_opt{'yAxis_TickUnitLen'}; 
			$show_v .= "$frame_opt{'yAxis_TickUnit'}"; 
			$grps{'yAxis_Tick'}->line(
				'id' => "yAxis_tick:$rI:$id:$i", 
				'x1' => $frame_lx,                             'y1' => -1 * ($i-1) / $frame_opt{'yAxis_BpPoint'}, 
				'x2' => $frame_lx-$frame_opt{'yAxis_TickLen'}, 'y2' => -1 * ($i-1) / $frame_opt{'yAxis_BpPoint'}, 
				'transform' => join('', "translate(0," , $frame_by - $sP[$rI], ")"), 
			); 
			my $xV = $frame_lx - $frame_opt{'yAxis_TickLen'}; 
			my $yV = -1 * ($i-1) / $frame_opt{'yAxis_BpPoint'}; 
			my $yTrans = $frame_by - $sP[$rI]; 
			$grps{'yAxis_Tick'}->text(
				'id' => "yAxis_tickTxt:$rI:$id:$i", 
				'x' => $xV , 'y' => $yV, 
				-cdata => $show_v, 
				'transform' => join('', "translate(0 $yTrans) rotate(90 $xV $yV)"), 
				'alignment-baseline'=> $frame_opt{'yAxis_TxtVertAln'}, 
				'stroke' => 'none', 
			); 
		}
	}# grps_yAxis_Tick : ticks
	my $chrID_xV = $frame_lx - $frame_opt{'yAxis_TickLen'} - $frame_opt{'yAxis_TickFontSize'} - $frame_opt{'yAxis_ChrIDFontSize'} * 0.3; 
	my $chrID_yV = $frame_by - ($sP[$rI] + $eP[$rI])/2 ; 
	$grps{'yAxis_Label'}->text(
		'x' => $chrID_xV, 'y' => $chrID_yV, 
		-cdata=>"$id", 
		'transform' => join('', "rotate(90 $chrID_xV $chrID_yV)"), 
		'alignment-baseline'=> $frame_opt{'yAxis_TxtVertAln'}, 
		'font-size' => $frame_opt{'yAxis_ChrIDFontSize'} , 
	); # For chromosome ID text 
}
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
# s04 - Draw syntenic blocks : 
##############################

for my $blkR (@tabInfo) {
	my ($blkID, $chrLoc_1, $chrLoc_2, $strand, $alnScore, $alnEvalue, $alnNumber, $geneLis_1, $geneLis_2, $kaLis, $ksLis, $kaksLis) = @$blkR; 
	if (defined $chrLisX{'repN'}{$chrLoc_1->[0]} and defined $chrLisY{'repN'}{$chrLoc_2->[0]}) {
		my @linkXY = &_seg2xy(
			'chrLisX' => \%chrLisX, 
			'chrLisY' => \%chrLisY, 
			'Ps' => [ $chrLoc_1->[0], $chrLoc_1->[1], $chrLoc_2->[0], $chrLoc_2->[1] ], 
			'Pe' => [ $chrLoc_1->[0], $chrLoc_1->[2], $chrLoc_2->[0], $chrLoc_2->[2] ], 
		); 
		$grps{'blk_line'}->line(
			'x' => , 'y' => , 
			''
		); 
	}
	if (defined $chrLisX{'repN'}{$chrLoc_2->[0]} and defined $chrLisY{'repN'}{$chrLoc_1->[0]}) {
	}
}

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

##############################
# sXX - Print out 
##############################
print {$outSvgFh} $svg->xmlify; 

################################################################################
#    Sub-routines 
################################################################################

################################################################################
#    Inner sub-routines 
################################################################################
sub _readInAddLine {
	my $fn = shift; 
	my %back; 
	defined $fn or return (\%back); 
	my $fh = &openFH( $fn, '<' ); 
	while (&wangLineC($fh)) {
		my @ta = &splitL("\t", $_); 
		push(@{$back{$ta[0]}}, $ta[1]); 
	}
	close($fh); 
	return(\%back); 
} # _readInAddLine() 
sub _shadow_setopt {
	for my $tk (qw/xAxis_shadow_setopt yAxis_shadow_setopt/) {
		defined $opts{$tk} or next; 
		for my $t1 (@{$opts{$tk}}) {
			for my $t2 ( split(/;/, $opts{$tk}{$t1}) ) {
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
					$ta[1] > $ta[2] and @ta[1,2] = @ta[2,1]; 
					push(@{$blk{$ta[0]}}, [@ta[1,2]]); 
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
	return; 
}# _blk_setopt() 


sub _img_setopt {
	&_setoptAR( \%img_opt, $opts{'img_setopt'}, 1 ); 
	&_setoptAR( \%img_opt, $default_string{'img_setopt'}, 0 ); 
	return; 
}# _img_setopt() 

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



