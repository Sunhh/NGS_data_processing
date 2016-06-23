package plotSunhh; 

use strict; 
use warnings; 
use Exporter qw(import);
our @EXPORT = qw(ins_calc);
our @EXPORT_OK = qw();

use LogInforSunhh; 
use mathSunhh; 

############################################################
#  Basic informations 
############################################################
my %html_to_rgb; # maroon => 'rgb(128,0,0)'
my %rgb_to_html; # 'rgb(139,0,0)' => 'dark_red'
{
	# Block for rgb table. 
	my $col_table = <<'CC'; 
maroon        #800000  (128,0,0)
dark_red      #8B0000  (139,0,0)
brown         #A52A2A  (165,42,42)
firebrick     #B22222  (178,34,34)
crimson       #DC143C  (220,20,60)
red           #FF0000  (255,0,0)
tomato        #FF6347  (255,99,71)
coral         #FF7F50  (255,127,80)
indian_red    #CD5C5C  (205,92,92)
light_coral   #F08080  (240,128,128)
dark_salmon   #E9967A  (233,150,122)
salmon        #FA8072  (250,128,114)
light_salmon  #FFA07A  (255,160,122)
orange_red    #FF4500  (255,69,0)
dark_orange   #FF8C00  (255,140,0)
orange        #FFA500  (255,165,0)
gold          #FFD700  (255,215,0)
dark_golden_rod  #B8860B  (184,134,11)
golden_rod       #DAA520  (218,165,32)
pale_golden_rod  #EEE8AA  (238,232,170)
dark_khaki       #BDB76B  (189,183,107)
khaki            #F0E68C  (240,230,140)
olive            #808000  (128,128,0)
yellow           #FFFF00  (255,255,0)
yellow_green     #9ACD32  (154,205,50)
dark_olive_green #556B2F  (85,107,47)
olive_drab       #6B8E23  (107,142,35)
lawn_green       #7CFC00  (124,252,0)
chart_reuse      #7FFF00  (127,255,0)
green_yellow     #ADFF2F  (173,255,47)
dark_green       #006400  (0,100,0)
green            #008000  (0,128,0)
forest_green     #228B22  (34,139,34)
lime             #00FF00  (0,255,0)
lime_green       #32CD32  (50,205,50)
light_green      #90EE90  (144,238,144)
pale_green       #98FB98  (152,251,152)
dark_sea_green   #8FBC8F  (143,188,143)
medium_spring_green  #00FA9A  (0,250,154)
spring_green         #00FF7F  (0,255,127)
sea_green            #2E8B57  (46,139,87)
medium_aqua_marine   #66CDAA  (102,205,170)
medium_sea_green     #3CB371  (60,179,113)
light_sea_green      #20B2AA  (32,178,170)
dark_slate_gray      #2F4F4F  (47,79,79)
teal                 #008080  (0,128,128)
dark_cyan            #008B8B  (0,139,139)
aqua                 #00FFFF  (0,255,255)
cyan                 #00FFFF  (0,255,255)
light_cyan           #E0FFFF  (224,255,255)
dark_turquoise       #00CED1  (0,206,209)
turquoise            #40E0D0  (64,224,208)
medium_turquoise  #48D1CC  (72,209,204)
pale_turquoise    #AFEEEE  (175,238,238)
aqua_marine       #7FFFD4  (127,255,212)
powder_blue       #B0E0E6  (176,224,230)
cadet_blue        #5F9EA0  (95,158,160)
steel_blue        #4682B4  (70,130,180)
corn_flower_blue  #6495ED  (100,149,237)
deep_sky_blue     #00BFFF  (0,191,255)
dodger_blue       #1E90FF  (30,144,255)
light_blue        #ADD8E6  (173,216,230)
sky_blue          #87CEEB  (135,206,235)
light_sky_blue    #87CEFA  (135,206,250)
midnight_blue     #191970  (25,25,112)
navy              #000080  (0,0,128)
dark_blue         #00008B  (0,0,139)
medium_blue       #0000CD  (0,0,205)
blue              #0000FF  (0,0,255)
royal_blue        #4169E1  (65,105,225)
blue_violet       #8A2BE2  (138,43,226)
indigo            #4B0082  (75,0,130)
dark_slate_blue   #483D8B  (72,61,139)
slate_blue        #6A5ACD  (106,90,205)
medium_slate_blue  #7B68EE  (123,104,238)
medium_purple  #9370DB  (147,112,219)
dark_magenta   #8B008B  (139,0,139)
dark_violet    #9400D3  (148,0,211)
dark_orchid    #9932CC  (153,50,204)
medium_orchid  #BA55D3  (186,85,211)
purple         #800080  (128,0,128)
thistle        #D8BFD8  (216,191,216)
plum           #DDA0DD  (221,160,221)
violet         #EE82EE  (238,130,238)
magenta  #FF00FF  (255,0,255)
fuchsia  #FF00FF  (255,0,255)
orchid   #DA70D6  (218,112,214)
medium_violet_red  #C71585  (199,21,133)
pale_violet_red  #DB7093  (219,112,147)
deep_pink        #FF1493  (255,20,147)
hot_pink    #FF69B4  (255,105,180)
light_pink  #FFB6C1  (255,182,193)
pink        #FFC0CB  (255,192,203)
antique_white  #FAEBD7  (250,235,215)
beige  #F5F5DC  (245,245,220)
bisque  #FFE4C4  (255,228,196)
blanched_almond  #FFEBCD  (255,235,205)
wheat  #F5DEB3  (245,222,179)
corn_silk  #FFF8DC  (255,248,220)
lemon_chiffon  #FFFACD  (255,250,205)
light_golden_rod_yellow  #FAFAD2  (250,250,210)
light_yellow  #FFFFE0  (255,255,224)
saddle_brown  #8B4513  (139,69,19)
sienna  #A0522D  (160,82,45)
chocolate  #D2691E  (210,105,30)
peru  #CD853F  (205,133,63)
sandy_brown  #F4A460  (244,164,96)
burly_wood  #DEB887  (222,184,135)
tan  #D2B48C  (210,180,140)
rosy_brown  #BC8F8F  (188,143,143)
moccasin  #FFE4B5  (255,228,181)
navajo_white  #FFDEAD  (255,222,173)
peach_puff  #FFDAB9  (255,218,185)
misty_rose  #FFE4E1  (255,228,225)
lavender_blush  #FFF0F5  (255,240,245)
linen  #FAF0E6  (250,240,230)
old_lace  #FDF5E6  (253,245,230)
papaya_whip  #FFEFD5  (255,239,213)
sea_shell  #FFF5EE  (255,245,238)
mint_cream  #F5FFFA  (245,255,250)
slate_gray  #708090  (112,128,144)
light_slate_gray  #778899  (119,136,153)
light_steel_blue  #B0C4DE  (176,196,222)
lavender  #E6E6FA  (230,230,250)
floral_white  #FFFAF0  (255,250,240)
alice_blue  #F0F8FF  (240,248,255)
ghost_white  #F8F8FF  (248,248,255)
honeydew  #F0FFF0  (240,255,240)
ivory  #FFFFF0  (255,255,240)
azure  #F0FFFF  (240,255,255)
snow  #FFFAFA  (255,250,250)
black  #000000  (0,0,0)
dim_gray  #696969  (105,105,105)
dim_grey  #696969  (105,105,105)
gray  #808080  (128,128,128)
grey  #808080  (128,128,128)
dark_gray  #A9A9A9  (169,169,169)
dark_grey  #A9A9A9  (169,169,169)
silver     #C0C0C0  (192,192,192)
light_gray  #D3D3D3  (211,211,211)
light_grey  #D3D3D3  (211,211,211)
gainsboro  #DCDCDC  (220,220,220)
white_smoke  #F5F5F5  (245,245,245)
white  #FFFFFF  (255,255,255)
CC

	for my $l1 ( split(/\n/, $col_table) ) {
		$l1 =~ m!^(\S+)\s+(#\S+)\s+\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)! or &stopErr("[Err] Failed to parse line [$l1]\n"); 
		my ($html_name, $hex_name, $rgb_1, $rgb_2, $rgb_3) = ($1, $2, $3,$4,$5); 
		$html_name = lc($html_name); 
		$hex_name  = lc($hex_name); 
		$html_to_rgb{$html_name} //= [$rgb_1,$rgb_2,$rgb_3]; 
		$rgb_to_html{"rgb($rgb_1,$rgb_2,$rgb_3)"} //= $html_name; 
	}
}# block for rgb_table. 


############################################################
#  Sub-routines
############################################################


=head1 _rgb_color( 'in_color'=>'html_name:black', 'out_fmt'=>'rgb' )

Return       : ($rgb_color)

  &_rgb_color( 'in_color'=>'html_name:black', 'out_fmt'=>'rgb' ) returns 'rgb(0,0,0)'
  &_rgb_color( 'in_color'=>'html_name:abcd' , 'out_fmt'=>'rgb' ) returns () # Error 
  &_rgb_color( 'in_color'=>'hex_name:#FFFFFF', 'out_fmt'=>'rgb' ) returns 'rgb(255,255,255)' # White 
	&_rgb_color( 'in_color'=>'html_name:Maroon', 'out_fmt'=>'hex' ) returns '#800000' 
	&_rgb_color( 'in_color'=>'rgb_name:rgb(240,128,128)', 'out_fmt'=>'html' ) returns 'light coral' 
	&_rgb_color( 'in_color'=>'html_name:light coral', 'out_fmt'=>'rgb' ) returns 'rgb(240,128,128)' 

=cut
sub _rgb_color {
	my %parm = &mathSunhh::_setHashFromArr(@_); 
	# $parm{'in_color'} = lc( $parm{'in_color'} ); 
	# $parm{'out_fmt'} = lc($parm{'out_fmt'}); 
	my @rgb; 
	if ( $parm{'in_color'} =~ m/^(?:html|css)_name:(.+)$/i ) {
		my $html_name = $1; 
		$html_name =~ s!\s+!_!g; $html_name = lc($html_name); 
		# defined $html_to_rgb{$html_name} or do { &tsmsg("[Err] Failed to find html_name [$html_name]\n"); return; }; 
		defined $html_to_rgb{$html_name} or return; 
		@rgb = @{$html_to_rgb{$html_name}}; 
	} elsif ( $parm{'in_color'} =~ m/^hex_name:#(\S{2})(\S{2})(\S{2})$/i ) {
		my ($h1, $h2, $h3) = ($1, $2, $3); 
		@rgb = &mathSunhh::_hexa_to_decimal( $h1, $h2, $h3 ); 
	} elsif ( $parm{'in_color'} =~ m/^rgb_name:rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i ) {
		@rgb = ($1, $2, $3); 
	} else {
		&stopErr("[Err] Failed to parse input in_color [$parm{'in_color'}]\n"); 
	}
	
	if ( $parm{'out_fmt'} =~ m/^rgb$/i ) {
		return("rgb($rgb[0],$rgb[1],$rgb[2])"); 
	} elsif ( $parm{'out_fmt'} =~ m/^hex$/i ) {
		return( '#' . join('', &mathSunhh::_decimal_to_hexa(@rgb)) );
	} elsif ( $parm{'out_fmt'} =~ m/^(html|css)$/i ) {
		my $tk = "rgb($rgb[0],$rgb[1],$rgb[2])"; 
		# defined $rgb_to_html{$tk} or do { &tsmsg("[Err] Failed to find html_name for [$tk]\n"); return; }; 
		defined $rgb_to_html{$tk} or return; 
		( my $back = $rgb_to_html{$tk} ) =~ s!_! !g; 
		return( $back ); 
	} else {
		&stopErr("[Err] Failed to parse out_fmt [$parm{'out_fmt'}], which should be one of 'rgb|hex|html';\n"); 
	}
	
	return; 
} # _rgb_color

=head1 cnvt_to_rgb ( $min, $max, $val, $col ) 

Function : Get a series of color for heat map. Invoke mathSunhh::cnvt_to_rgb(); 

Return : ( "rgb(R,G,B)" )

Required : 
  $col = [ [r1,g1,b1], [r2,g2,b2], [r3,g3,b3], ... ]
  $min , $max : value range; 
  $val : current value. 

Reference: http://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map

=cut
sub cnvt_to_rgb {
	return( &mathSunhh::cnvt_to_rgb(@_) ); 
#	my ( $min, $max, $val, $col ) = @_; 
#	my $col_i = $#$col;
#	my $v = ($val-$min)/($max-$min) * $col_i;
#	my $i1 = &min( int($v), $col_i );
#	my $i2 = &min( int($v)+1, $col_i);
#	my $f = $v - $i1;
#	my @back;
#	for (my $j=0; $j<3; $j++) {
#		$back[$j] = int( $col->[$i1][$j] + $f * ($col->[$i2][$j]-$col->[$i1][$j]) );
#	}
#	return( sprintf("rgb(%.0f,%.0f,%.0f)", $back[0], $back[1], $back[2])  );
}# cnvt_to_rgb () 



1;