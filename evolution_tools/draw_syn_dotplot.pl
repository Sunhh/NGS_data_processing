#!/usr/bin/perl
# http://search.cpan.org/~darnold/DBD-Chart-0.82/Chart/Plot.pm#Establish_data_points:_setPoints() 
use strict; 
use warnings; 
use LogInforSunhh; 
use fileSunhh; 
use DBD::Chart::Plot;
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	# Global 
	"img_width:i", 
	"img_height:i", 

	"bp_per_point:i", 

	"out_fmt:s", # png 
	"out_name:s", # out.$opts{'out_fmt'}

	"img_setopt:s", # 'title=;horizMargin=50;vertMargin=50;xAxisLabel=Genome 1;yAxisLabel=Genome 2'
	"no_scf_border!", # Do not add scaffolds' border in figure. 
	"add_x_border:s", "add_y_border:s", # Add border lines for x/y axises 

	# Scaffolds list 
	"scfLis_x:s", 
	"scfLis_y:s", 
	# gene_position list 
	"in_gff:s", 
	# mcscan_out.collinearity 
	"in_aln:s", 
); 

$opts{'bp_per_point'} //= 100; # bps number each point in figure standing for. 
$opts{'out_fmt'}  //= 'png'; 
$opts{'out_name'} //= "out.$opts{'out_fmt'}"; 
$opts{'img_setopt'} //= 'title=;horizMargin=50;vertMargin=50;xAxisLabel=Genome 1;yAxisLabel=Genome 2'; 

my $help_txt = <<HH; 

perl $0    -scfLis_x   scflist.x     -scfLis_y   scflist.y    -in_gff  mcs.gff    -in_aln  mcs.collinearity

-help 

-img_width       [50 + max_len/bp_per_point]
-img_height      [50 + max_len/bp_per_point]
-bp_per_point    [$opts{'bp_per_point'}]

-img_setopt      [$opts{'img_setopt'}]

-no_scf_border   [Boolean]
-add_x_border    [filename] Add border at the end of ScfID . Format: ScfID|position \\n
-add_y_border    [filename] Please DON'T use number as ScfID!!!!

-out_fmt         [$opts{'out_fmt'}]
-out_name        [$opts{'out_name'}]

#### Format of scfLis_? :   Scf_ID  \\t  Scf_Len  
#### Format of in_gff   :   Scf_ID  \\t  Gene_ID  \\t  Gene_Start  \\t  Gene_End (unstranded)
#### Format of in_aln   :   Output of mcscanx 
#### The Scf_IDs in scfList_? must be same to those in in_gff, but not necessarily same to those in in_aln file. 
#### The gene_IDs in in_aln file must exist in in_gff files. 


HH

$opts{'help'} and &LogInforSunhh::usage($help_txt); 

######################################################################
#    Reading data 
######################################################################
my ($chr2gen, $gen2loc) = &_readInGff($opts{'in_gff'}) ; 
my $alnInfo = &_readInAln( $opts{'in_aln'} ); 
defined $alnInfo->[0]{'info'} or shift(@{$alnInfo}); 
my ($scfInfo_x) = &_readInScfLis( $opts{'scfLis_x'} ); 
my ($scfInfo_y) = &_readInScfLis( $opts{'scfLis_y'} ); 
my $ofh = &openFH($opts{'out_name'}, '>'); 
my %border_list; # {x|y} => {border_IDs => 1}; 
for (@{$scfInfo_x->{'arr'}}) { $border_list{'x'}{$_->[0]} = $_->[0]; }
for (@{$scfInfo_y->{'arr'}}) { $border_list{'y'}{$_->[0]} = $_->[0]; }
$opts{'no_scf_border'} and do { $border_list{'x'} = {}; $border_list{'y'} = {}; }; 
if (defined $opts{'add_x_border'}) {
	for my $t1 (@{ &_readInBorder($opts{add_x_border}) }) { $border_list{'x'}{$t1} = $t1; }
}
if (defined $opts{'add_y_border'}) {
	for my $t1 (@{ &_readInBorder($opts{add_y_border}) }) { $border_list{'y'}{$t1} = $t1; }
}
keys %{$border_list{'x'}} > 0 and &_borderAsNumber( $border_list{'x'}, $scfInfo_x ); 
keys %{$border_list{'y'}} > 0 and &_borderAsNumber( $border_list{'y'}, $scfInfo_y ); 

######################################################################
#    Compute img parameters 
######################################################################
my %img_para; 
$img_para{'img_width'}  //= $opts{'img_width'}  ; 
$img_para{'img_height'} //= $opts{'img_height'} ; 
$img_para{'max_x'} = $scfInfo_x->{'arr'}[-1][2] + $scfInfo_x->{'arr'}[-1][1]; 
$img_para{'max_y'} = $scfInfo_y->{'arr'}[-1][2] + $scfInfo_y->{'arr'}[-1][1]; 
unless (defined $img_para{'img_width'}) {
	$img_para{'img_width'}  = 50 + $img_para{'max_x'} / $opts{'bp_per_point'}; 
&tsmsg("width = $img_para{'img_width'}\n"); 
}
unless (defined $img_para{'img_height'}) {
	$img_para{'img_height'} = 50 + $img_para{'max_y'} / $opts{'bp_per_point'}; 
}

my $img = DBD::Chart::Plot->new( $img_para{'img_width'} , $img_para{'img_height'} ); 
&set_img_opt( $img, $opts{'img_setopt'} ); 
&set_border( $img, $border_list{'x'}, $img_para{'max_y'}, 'x' ); 
&set_border( $img, $border_list{'y'}, $img_para{'max_x'}, 'y' ); 

## Setup the whole size of x-y plot 
# $img->setPoints( [ 1, $img_para{'max_x'} ], [ 1, $img_para{'max_y'} ], 'black noline nopoints dot' ); 
$img->setPoints( [ $img_para{'max_x'}, $img_para{'max_x'} ], [ 1,                  $img_para{'max_y'} ], 'black line nopoints dot' ); 
$img->setPoints( [ 1,                  $img_para{'max_x'} ], [ $img_para{'max_y'}, $img_para{'max_y'} ], 'black line nopoints dot' ); 

## For the points in each alignment , set the points-in-line and add them to $img ; 
for (my $i=0; $i<@{$alnInfo}; $i++) {
	my (@gen1, @gen2, @ka, @ks, @w); 
	for my $ar1 (@{$alnInfo->[$i]{'pair'}}) {
		push(@gen1, $ar1->[0]); 
		push(@gen2, $ar1->[1]); 
		push(@ka, $ar1->[3]); 
		push(@ks, $ar1->[4]); 
		push(@w,  $ar1->[5]); 
	}
	my ($alnID, $score, $eval, $npair, $chr1, $chr2, $str) = @{$alnInfo->[$i]{'info'}}; 
	$str eq 'plus'  and $str = '+'; 
	$str eq 'minus' and $str = '-'; 

	# Use the chrID in in_gff file. 
	$chr1 = $gen2loc->{ $gen1[0] }[0]; 
	$chr2 = $gen2loc->{ $gen2[0] }[0]; 
	
	my $color = 'blue'; 
	$str eq '-' and $color = 'purple'; 

	my (@xdata, @ydata); 
	if ( defined $scfInfo_x->{'hash'}{$chr1} and defined $scfInfo_y->{'hash'}{$chr2} ) {
		for ( my $i=0; $i<@gen1; $i++ ) {
			for my $j ( @{ $scfInfo_x->{'hash'}{$chr1} } ) {
				# For X 
				my $x_base = $scfInfo_x->{'arr'}[$j][2]; 
				my $x_v_s = $x_base + $gen2loc->{ $gen1[$i] }[1]; 
				my $x_v_e = $x_base + $gen2loc->{ $gen1[$i] }[2]; 
				for my $k ( @{ $scfInfo_y->{'hash'}{$chr2} } ) {
					# For Y 
					my $y_base = $scfInfo_y->{'arr'}[$k][2]; 
					my $y_v_s = $y_base + $gen2loc->{ $gen2[$i] }[1]; 
					my $y_v_e = $y_base + $gen2loc->{ $gen2[$i] }[2]; 
					push(@xdata, $x_v_s); 
					push(@ydata, $y_v_s); 
				}
			}
		}
		$img->setPoints(\@xdata, \@ydata, $color . ' line points dot'); 
	}
	@xdata=(); @ydata=(); 
	if ( defined $scfInfo_x->{'hash'}{$chr2} and defined $scfInfo_y->{'hash'}{$chr1} ) {
		for ( my $i=0; $i<@gen2; $i++ ) {
			for my $j ( @{ $scfInfo_x->{'hash'}{$chr2} } ) {
				# For X 
				my $x_base = $scfInfo_x->{'arr'}[$j][2]; 
				my $x_v_s = $x_base + $gen2loc->{ $gen2[$i] }[1]; 
				my $x_v_e = $x_base + $gen2loc->{ $gen2[$i] }[2]; 
				for my $k ( @{ $scfInfo_y->{'hash'}{$chr1} } ) {
					# For Y 
					my $y_base = $scfInfo_y->{'arr'}[$k][2]; 
					my $y_v_s = $y_base + $gen2loc->{ $gen1[$i] }[1]; 
					my $y_v_e = $y_base + $gen2loc->{ $gen1[$i] }[2]; 
					push(@xdata, $x_v_s); 
					push(@ydata, $y_v_s); 
				}
			}
		}
		$img->setPoints(\@xdata, \@ydata, $color . ' line points dot'); 
	} else {
		# Next; 
	}
}

print {$ofh} $img->plot($opts{'out_fmt'}); 
close($ofh); 


######################################################################
#     Sub-routines 
######################################################################

# Return : {} 
#   'arr' => ([ [scf1_ID, scf1_Len, cum_len_prevE], [scf2_ID, scf2_Len, cum_len_prevE], ... ]); 
#   'hash' => { scf_ID => [idx_in_arr, ...] }
sub _readInScfLis {
	# $_[0] : input file name; Format: Scf_ID \\t  Scf_Len  
	my $inFh = &openFH($_[0], '<'); 
	my %back; 
	my $cum_len = 0; 
	while (<$inFh>) {
		chomp; 
		m/^\s*(#|$)/ and next; 
		my @ta = split(/\t/, $_); 
		( defined $ta[1] and $ta[1] > 0 ) or do { &tsmsg("[Wrn] Skip scaffold [$ta[0]] with bad length [$ta[1]].\n"); next; }; 
		defined $back{'hash'}{$ta[0]} and &tsmsg("[Wrn] Repeated scaffold ID [$ta[0]]\n"); 
		push(@{$back{'arr'}}, [ $ta[0], $ta[1], $cum_len ]); 
		push(@{ $back{'hash'}{$ta[0]} }, $#{ $back{'arr'} } ); 
		$cum_len += $ta[1]; 
	}
	close($inFh); 
	return(\%back); 
}# _readInScfLis() 


# Return : (\@alnInfo)
#  @alnInfo : 
#    [idx_num]{'info'} => [alnID, score, evalue, Num_pairs, chrID_1, chrID_2, strand(plus|minus)]
#    [idx_num]{'pair'} => [ [genID_1, genID_2, pair_evalue, Ka, Ks, Ka/Ks], [], ... ]
#    [idx_num]{'text'} => $text_of_current_block
sub _readInAln{
	my $inFh = &openFH(shift, '<'); 
	my @alnInfo; 
	my $aln_idx = 0; 
	while (<$inFh>) {
		chomp; 
		if ( m!^\s*$|^####|^# \S+:\s+\S+$|# (MAX GAPS|Number of collinear genes|Number of all genes):! ) {
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
			next; 
		}
		if ( m!^## Alignment ! ) {
			# ## Alignment 0: score=7137.0 e_value=0 N=147 Ma1&Ma1 plus
			m!## Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) ([^\&\s]+)\&([^\&\s]+) (plus|minus|X+)! or die "ALN:$_\n"; 
			my ($alnID, $score, $eval, $n, $chr1, $chr2, $str)
			=  ($1,     $2,     $3,    $4, $5,    $6,    $7); 
			$aln_idx++; 
			$alnInfo[$aln_idx]{'info'} = [$alnID, $score, $eval, $n, $chr1, $chr2, $str]; 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} elsif ( m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+))?)?$! ) { 
			# #  0-  0:        Cma_000007      Cma_000973        2e-57
			$aln_idx > -1 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tka, $tks, $tw) 
			= 
			   ($1,     $2,        $3,    $4,    $5,    $6,   $7,   $8); 
			$tka //= ''; $tks //= ''; $tw //= ''; 
			$alnID == $alnInfo[$aln_idx]{'info'}[0] or &stopErr("[Err] line_alnID=$alnID not fitting upper level (alnID=$alnInfo[$aln_idx]{'info'}[0]).\n"); 
			push(@{$alnInfo[$aln_idx]{'pair'}}, [$gid1, $gid2, $eval, $tka, $tks, $tw]); 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} elsif ( m!^\s*(\d+)\-\s*(\d+):\t(\S+)\t(\S+)\s*(\S+)(?:\t(\S+)\t(\S+)(?:\t(\S+)\t(\S+))?)?$! ) { 
			# #  0-  0:  Cma_000007      Cma_000973        2e-57 0.5258  0.0154  0.5387  0.0148
			# This is to fit result from haibao tang's python ks calculation. 
			$aln_idx > -1 or &stopErr("Too early to line: $_\n"); 
			my ($alnID, $alnID_id, $gid1, $gid2, $eval, $tks, $tka, $t_ngKs, $t_ngKa) 
			= 
			(   $1,     $2,        $3,    $4,    $5,    $6,   $7,   $8,      $9); 
			my $tw; 
			if ( defined $t_ngKs ) {
				$tks = $t_ngKs; 
				$tka = $t_ngKa; 
			}
			$tks < 0 and $tks = 'nan'; 
			if ( defined $tks and $tks ne 'nan') {
				$tw = ( $tks > 0 ) ? $tka/$tks : 'nan'; 
			}
			$tks eq 'nan' and $tw = 'nan'; 
			$tka //= ''; $tks //= ''; $tw //= ''; 
			$alnID == $alnInfo[$aln_idx]{'info'}[0] or &stopErr("[Err] line_alnID=$alnID not fitting upper level (alnID=$alnInfo[$aln_idx]{'info'}[0]).\n"); 
			push(@{$alnInfo[$aln_idx]{'pair'}}, [$gid1, $gid2, $eval, $tka, $tks, $tw]); 
			$alnInfo[$aln_idx]{'text'} .= "$_\n"; 
		} else {
			&stopErr("[Err] Unable to parse line: $_\n"); 
		}
	}
	close($inFh); 
	return(\@alnInfo); 
}# _readInAln

# Return : (\%chr2gen, \%gen2loc)
#  %chr2gen : {chromID} => [ [genID, chrS, chrE], [], ... ]
#  %gen2loc : {genID} => [chrID, chrS, chrE]
sub _readInGff {
	my $inFh = &openFH(shift, '<'); 
	my (%chr2gen, %gen2loc); 
	while (<$inFh>) {
		chomp; 
		my ( $chrID, $genID, $chrS, $chrE ) = split(/\t/, $_); 
		push( @{$chr2gen{$chrID}}, [$genID, $chrS, $chrE] ); 
		defined $gen2loc{$genID} and &tsmsg("[Msg] genID [$genID] repeats, and loc infor masked by latter one\n"); 
		$gen2loc{$genID} = [ $chrID, $chrS, $chrE ]; 
	}
	close($inFh); 
	return (\%chr2gen, \%gen2loc); 
}# _readInGff() 

# Function : _borderAsNumber( $border_list{'x'} , $scfInfo_x )
#  Convert values in $border_list{'x'}{$border} to numbers according to $scfInfo_x 
#  Please don't use number as scfID in $scfInfo_x !!!! 
sub _borderAsNumber {
	# $_[0] : $border_list{'x'} 
	# $_[1] : $scfInfo_x 
	my @bad_border; 
	for my $k ( keys %{ $_[0] } ) {
		if ( defined $_[1]->{'hash'}{$k} ) {
			my $i = $_[1]->{'hash'}{$k}[0]; 
			$_[0]->{ $k } = $_[1]->{'arr'}[$i][2]+$_[1]->{'arr'}[$i][1]; 
		} elsif ( $k =~ m/^\-?[\d\.]+$/i ) {
			$_[0]->{ $k } = $k; 
			; 
		} else {
			&tsmsg("[Wrn] Remove bad assignment [$k] for border file\n"); 
			push(@bad_border, $k); 
		} 
	}
	for my $k (@bad_border) {
		delete $_[0]->{$k}; 
	}
	return undef(); 
}# _borderAsNumber() 

# Function : _readInBorder( add_x_border.filename ) 
#  Format of add_x_border.filename : number \\n scfID \\n ... 
# Return   : [ the_first_col ]
sub _readInBorder {
	# $_[0] : add_x_border.filename
	my $fh = &openFH($_[0], '<'); 
	my %back; 
	while (<$fh>) {
		chomp; m/^\s*(#|$)/ and next; 
		my @ta = split(/\t/, $_); 
		$back{$ta[0]} //= $.; 
	}
	close($fh); 
	return ([sort { $back{$a} <=> $back{$b} } keys %back]); 
}# _readInBorder() 

sub set_img_opt {
	# $_[0] : $img object 
	# $_[1] : $opts{'img_setopt'}
	my %para_1 = (
		'horizMargin'  =>  75, 
		'vertMargin'   =>  50, 
		'title'        =>  '', 
		'xAxisLabel'   =>  'Genome 1', 
		'yAxisLabel'   =>  'Genome 2', 
	); 
	for my $t1 (split(/;/, $_[1])) {
		$t1 =~ s/(^\s+|\s+$)//g; 
		$t1 =~ s/^([^=]+)=// or do { &tsmsg("[Wrn] Skip bad settings [$t1] in -img_setopt\n"); next; }; 
		my $tag = $1; $tag =~ s/\s+$//; 
		$t1 =~ s/^\s+//; 
		$para_1{$tag} = $t1; 
	}
	$_[0]->setOptions( %para_1 ); 
	return undef(); 
}# set_img_opt 

sub set_border {
	# $_[0] : $img 
	# $_[1] : $border_list{'x'} as numeric 
	# $_[2] : $max_of_the_other_axis 
	# $_[3] : x | y 
	$_[3] //= 'x'; 
	for my $b1 (keys %{$_[1]}) {
		my @data1; 
		my @data2; 
		push(@data1, $_[1]->{$b1}, $_[1]->{$b1}); 
		push(@data2, 1,            $_[2]); 
		if ( $_[3] =~ m/^x$/i ) {
			$_[0]->setPoints( \@data1, \@data2, 'black line nopoints' ); 
		} elsif ( $_[3] =~ m/^y$/i ) {
			$_[0]->setPoints( \@data2, \@data1, 'black line nopoints' ); 
		} else {
			&stopErr("[Err] Bad axis [$_[3]]\n"); 
		}
	}
	return undef(); 
}# set_border() 

