#!/usr/bin/perl -w 
use strict; 
use LogInforSunhh; 
use fileSunhh; 
use mathSunhh; 
my $ms = mathSunhh->new(); 
use SVG; 
use Getopt::Long; 
my %opts; 
GetOptions(\%opts, 
	"help!", 
	"bam:s",   # required
	"scfID:s", # required 
	"outSvg:s", # out.svg 
	"scfS:i",  # default 1
	"scfE:i",  # required if -scfS 
	"bp_per_point:f", # default 100
	"width_per_line:i", # default 1000 
	"bp_overlap:i", # default 30000 
	"tickStep:f", # default 10000 
	"tickUnit:s", # default k 
	"gapLis:s",   # N gap list. 
); 

sub usage {
	print <<HH; 
################################################################################
# perl $0 -bam in.srt.bam -scfID scfID 
# 
# -bam            [] Required. input sorted bam file. 
# -scfID          [] Required. 
# -outSvg         [] Output svg text to this file instead of stdout if given. 
# -scfS           [1] 
# -scfE           [] End of scaffold
# -bp_per_point   [100] bps number each point in figure standing for. 
# -width_per_line [1000] Figure width per line. 
# -bp_overlap     [20e3]
# -tickStep       [10000]
# -tickUnit       [k] k/m/''
#
# -gapLis         [] NSP306_Pla01s04GC_Gt5h.scf.fa.Nlis 
################################################################################
HH
	exit 1; 
}
$opts{'help'} and &usage(); 
( defined $opts{'bam'} and defined $opts{'scfID'} ) or &usage(); 

$opts{'scfS'} //= 1; 
$opts{'bp_per_point'} //= 100; 
$opts{'width_per_line'} //= 1000; 
$opts{'bp_overlap'} //= 20e3; 
$opts{'tickStep'} //= 10000; 
$opts{'tickUnit'} //= 'k'; 
my $unit_len = 1000; 
$opts{'tickUnit'} eq '' and $unit_len = 1; 
$opts{'tickUnit'} eq 'k' and $unit_len = 1e3; 
$opts{'tickUnit'} eq 'm' and $unit_len = 1e6; 

my $outFh = \*STDOUT; 
defined $opts{'outSvg'} and $outFh = &openFH($opts{'outSvg'}, '>'); 


my @add_blks; 
if (defined $opts{'gapLis'}) {
	open G,'<',"$opts{'gapLis'}" or die; 
	while (<G>) {
		chomp; 
		my @ta = split(/\t/, $_); 
		my ($id, $s, $e) = @ta[0, 2, 3]; 
		$id eq $opts{'scfID'} or next; 
		$s > $e and ($s, $e) = ($e, $s); 
		push(@add_blks, [$s, $e]); 
	}
	close G; 
}


# Read infor from bam file. 
-e "$opts{'bam'}.bai" or &exeCmd_1cmd("samtools index $opts{'bam'}"); 
my @pair_se; 
{
my $loc = ( defined $opts{'scfE'} ) ? "$opts{'scfID'}:$opts{'scfS'}-$opts{'scfE'}" : "$opts{'scfID'}" ; 
open F, '-|',"samtools view $opts{'bam'} $loc | sam_filter.pl -h2diff_F " or die; 
while (<F>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	my ($id1, $pos1, $id2, $pos2, $ins_len) = @ta[2,3,6,7,8]; 
	$id2 eq '=' or next; 
	my $pos3 = $pos1+$ins_len-1; 
	$pos1 <= $pos3 or next; # Skip PE pairs. 
	$pos1 >= $opts{'scfS'} or next; 
	defined $opts{'scfE'} and $pos3 > $opts{'scfE'} and next; 
	push(@pair_se, [$pos1, $pos3]); 
}
close F; 
@pair_se = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @pair_se; 
$opts{'scfE'} //= $pair_se[-1][1]; 
}

# For svg: 
# Basic parameters 
my $base_y = 10; 
my $base_x = 100; 
my $step_y = 100; 
my $rdHeight_y = $step_y * 2/3; 
my $bp_per_line = int( $opts{'width_per_line'} * $opts{'bp_per_point'} ); 
&tsmsg("[Rec] bp_per_line=$bp_per_line\n"); 

my $width; 
my $height; 
my %bp_line_wind; 
$width  = $base_x * 2 + $opts{'width_per_line'}; 
my $max_idx_dep = ($opts{'scfE'}-$opts{'scfS'}+1)/($bp_per_line-$opts{'bp_overlap'}); 
$max_idx_dep == int($max_idx_dep) or $max_idx_dep = int($max_idx_dep+1); 
my $base_y_h = ( $base_y >= 50 ) ? $base_y : 50 ; 
$height = $base_y + $base_y_h + $step_y * $max_idx_dep; 
%bp_line_wind = %{ 
 $ms->setup_windows( 
   'ttl_start' => $opts{'scfS'}, 
   'ttl_end'   => $opts{'scfE'}, 
   'wind_size' => $bp_per_line, 
   'wind_step' => $bp_per_line - $opts{'bp_overlap'}, 
   'minRatio'  => ($opts{'bp_overlap'}+1)/$bp_per_line, 
 ) 
}; 
my %si_to_idx; # {si} => index_of_line. 
for ( my $i=0; $i<@{$bp_line_wind{'info'}{'windSloci'}}; $i++ ) {
	$si_to_idx{ $bp_line_wind{'info'}{'windSloci'}[$i] } = $i+1; 
}


# SVG objects. 
my $svg = SVG->new('width'=>$width, 'height'=>$height); 
my %grps; 
$grps{'backbone'} = $svg->group(
	'id'                => "backbone_$opts{'scfID'}", 
	'stroke-width'      => 1, 
	'stroke'            => 'black', 
	'opacity'           => 1, 
	'text-anchor'       => 'middle', 
	'font-weight'       => 'normal', 
	'font-size'         => '10', 
	'font-family'       => 'ArialNarrow'
); # This is used to draw backbone of scaffold. 
$grps{'add_blks'} = $svg->group(
	'id'                => "add_blks_$opts{'scfID'}", 
	'stroke-width'      => 1, 
	'stroke'            => 'red', 
	'fill'              => 'red', 
	'opacity'           => 1, 
	'text-anchor'       => 'middle', 
	'font-weight'       => 'normal', 
	'font-size'         => '10', 
	'font-family'       => 'ArialNarrow'
); # This is used to draw add_blocks of scaffold. 
$grps{'readpair'} = $svg->group(
	'id'                => "read_pairs", 
	'stroke-width'      => 0.5, 
	'stroke'            => 'orange', 
	'opacity'           => 0.7, 
	'fill'              => 'transparent', 
	'text-anchor'       => 'middle', 
	'font-weight'       => 'normal', 
	'font-size'         => '10', 
	'font-family'       => 'ArialNarrow'
); # Draw read pairs. 

# Raw scaffoldID ; 
$grps{'backbone'}->text(
 'x'       => $base_x/2, 
 'y'       => $base_y/2, 
 -cdata    => "scfID=[$opts{'scfID'}] $opts{'scfS'}-$opts{'scfE'}", 
 'font-weight' => "bold", 
); 

# Draw ticks. 
for (my $i=0; $i<=$opts{'scfE'}; $i+=$opts{'tickStep'}) {
	$i >= $opts{'scfS'} or next; 
	my @si = @{ 
	 $ms->map_windows( 
	   'position'  => $i, 
	   'wind_hash' => \%bp_line_wind, 
	 ) 
	}; 
	my $show_v = $i/$unit_len; 
	$show_v .= $opts{'tickUnit'}; 
	for my $tsi (@si) {
		my $idx_dep = $si_to_idx{$tsi}; 
		$grps{'backbone'}->line(
		 'x1' => $base_x+($i-$bp_line_wind{'loci'}{$tsi}[0]+1)/$opts{'bp_per_point'}, 
		 'y1' => $base_y+$idx_dep*$step_y, 
		 'x2' => $base_x+($i-$bp_line_wind{'loci'}{$tsi}[0]+1)/$opts{'bp_per_point'}, 
		 'y2' => $base_y+$idx_dep*$step_y+5, 
		); 
		$grps{'backbone'}->text(
		 'x'       => $base_x+($i-$bp_line_wind{'loci'}{$tsi}[0]+1)/$opts{'bp_per_point'}, 
		 'y'       => $base_y+$idx_dep*$step_y+15, 
		 -cdata    => "$show_v", 
		 'font-weight' => 'lighter', 
		); 
	}
}

# Draw backbone lines. 
for ( my $i=1; $i<=$opts{'scfE'}; $i+=($bp_per_line - $opts{'bp_overlap'}) ) {
	my $lineE = $i+$bp_per_line-1; 
	my $lineS = $i; 
	$lineE < $opts{'scfS'} and next; 
	$lineE > $opts{'scfE'} and $lineE = $opts{'scfE'}; 
	$lineS < $opts{'scfS'} and $lineS = $opts{'scfS'}; 
	my @si = @{ 
	 $ms->map_windows( 
	   'position'  => $lineS, 
	   'wind_hash' => \%bp_line_wind, 
	 ) 
	}; 
	for my $tsi (@si) {
		my $idx_dep = $si_to_idx{$tsi}; 
		$grps{'backbone'}->line(
		 'x1' => $base_x, 
		 'y1' => $base_y+$idx_dep*$step_y, 
		 'x2' => $base_x+($lineE-$bp_line_wind{'loci'}{$tsi}[0]+1)/$opts{'bp_per_point'}, 
		 'y2' => $base_y+$idx_dep*$step_y, 
		); 
	}
}

# Draw add_blks lines. 
my %used_id; 
for my $tr (@add_blks) {
	my ($s, $e) = @$tr; 
	$e < $opts{'scfS'} and next; 
	$s > $opts{'scfE'} and next; 
	$s >= $opts{'scfS'} or $s = $opts{'scfS'}; 
	$e <= $opts{'scfE'} or $e = $opts{'scfE'}; 
	for ( my $i=1; $i<=$e; $i+=($bp_per_line - $opts{'bp_overlap'}) ) {
		my $lineE = $i+$bp_per_line-1; 
		my $lineS = $i; 
		$lineE < $s and next; 
		$lineE > $e and $lineE = $e; 
		$lineS < $s and $lineS = $s; 
		my @si = @{ 
		 $ms->map_windows( 
		   'position'  => $lineS, 
		   'wind_hash' => \%bp_line_wind, 
		 ) 
		}; 
		for my $tsi (@si) {
			my $idx_dep  = $si_to_idx{$tsi}; 
			my $cur_s    = $lineS-$bp_line_wind{'loci'}{$tsi}[0]+1; 
			my $cur_e    = $lineE-$bp_line_wind{'loci'}{$tsi}[0]+1; 
			my $cur_x_s  = $base_x + $cur_s/$opts{'bp_per_point'}; 
			my $cur_x_e  = $base_x + $cur_e/$opts{'bp_per_point'}; 
			my $tk = "$lineS - $lineE"; 
			my $n = 0; 
			while (defined $used_id{$tk}) {
				$tk = "$lineS - $lineE : $n"; 
				$n++; 
				$n > 10000 and die "Problem tk=$tk\n"; 
			}
			$grps{'add_blks'}->rectangle(
			 'x'      => $cur_x_s, 
			 'y'      => $base_y+$idx_dep*$step_y-8, 
			 'width'  => $cur_x_e-$cur_x_s, 
			 'height' => 5, 
			 'id'     => "$tk", 
			); 
			$used_id{$tk} = 1; 
		}
	}
}


# Draw read pairs. 
for my $tr (@pair_se) {
	my ($s, $e) = @$tr; 
	my @si_s = @{
	 $ms->map_windows(
	   'position'  => $s, 
	   'wind_hash' => \%bp_line_wind, 
	 )
	}; 
	my @si_e = @{
	 $ms->map_windows(
	   'position'  => $e, 
	   'wind_hash' => \%bp_line_wind, 
	 )
	}; 
	for my $tsi_s (@si_s) {
		my $idx_dep = $si_to_idx{$tsi_s}; 
		my $cur_s = $s-$bp_line_wind{'loci'}{$tsi_s}[0]+1; 
		my $cur_x_s  = $base_x + $cur_s/$opts{'bp_per_point'}; 
		for my $tsi_e (@si_e) {
			$tsi_s == $tsi_e or next; 
			my $cur_e = $e-$bp_line_wind{'loci'}{$tsi_e}[0]+1; 
			my $cur_x_e  = $base_x + $cur_e/$opts{'bp_per_point'}; 
			my $cur_y_se = $base_y+$idx_dep*$step_y; 
			my @xx = ( $cur_x_s,  ($cur_x_s+$cur_x_e)/2, $cur_x_e ); 
			my @yy = ( $cur_y_se-10, $cur_y_se-$rdHeight_y, $cur_y_se-10 ); 
			my $points = $grps{'readpair'}->get_path(
				'x'   => \@xx, 
				'y'   => \@yy, 
				-relative=>1, 
				-type=>'polyline',
				-closed=>0
			); 
			$grps{'readpair'}->polyline(%$points); 
		}
	}
}

print {$outFh} $svg->xmlify; 



