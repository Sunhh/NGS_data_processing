#!/usr/bin/perl -w 
use strict; 
use SVG; 
!@ARGV and die "perl $0 RIL_tag input.cov.1k\n"; 

my $tag = shift; 
my $covF = "A11_wind1k.cov.1k"; 
$covF = shift; 
(my $wind_size = $covF) =~ s/^.*\.(\d+)k$/$1/i; 

open CF,'<',"$covF" or die; 
my %hash; 
my %max_len; 
while (<CF>) {
  chomp; 
  my @ta = split(/\t/, $_); 
  # $ta[1] =~ /^\d+$/ or next; 
  push(@{$hash{$ta[0]}}, [@ta]); 
  if (defined $max_len{$ta[0]}) {
    $max_len{$ta[0]} < $ta[1] and $max_len{$ta[0]} = $ta[1]; 
  }else{
    $max_len{$ta[0]} = $ta[1]; 
  }
}
close CF; 

my $width = 3500; 
my $height = 600; 
my $y_base = 500; 
$height = 400; 
$y_base = 350; 
my $x_base = 100; 
my $colWid = 0.1; 
$colWid = 0.5; 
my $font_size = 15; 
my $xscale = 35000000/3500; 
my $yscale = 0.1; #  Number per point. 
$yscale = 0.5; # For 10kb window. 
if (defined $wind_size and $wind_size > 0) {
  $yscale = 0.5*($wind_size/10);
}

my %grpCol = qw(
homoPI 3 
hete 4 
homo97 5 
WinS 1 
); 

for my $chrn (sort keys %hash) {
  $chrn eq 'Chr' and next; 
 #  $chrn eq 'Chr2' or next; 
  my $svg = SVG->new(width=>$x_base+$width,height=>$height);
  my %grps; 
  $grps{homoPI} = $svg->group(
    id => 'homoPI', 
    'stroke' => 'red', 'fill' => 'red', 
    'stroke-width' => $colWid, 
  ); 
  $grps{hete} = $svg->group(
    id => 'hete', 
    'stroke' => 'green', 'fill' => 'green', 
    'stroke-width' => $colWid, 
  ); 
  $grps{homo97} = $svg->group(
    id => 'homo97', 
    'stroke' => 'blue', 'fill' => 'blue', 
    'stroke-width' => $colWid, 
  ); 
  $grps{Xbase} = $svg->group(
    id => 'XbaseLine', 
    'stroke-width' => 1, 
    'stroke' => 'black', 
    'text-anchor' => 'middle', 
    'font-weight' => "bold", 
    'font-size'   => $font_size, 
    'font-family' => "ArialNarrow", 
  ); 
  $grps{Ybase} = $svg->group(
    id => 'YbaseLine', 
    'stroke-width' => 1, 
    'stroke' => 'black', 
    'text-anchor' => 'end', 
    'font-weight' => "bold", 
    'font-size'   => $font_size, 
    'font-family' => "ArialNarrow", 
  ); 
  
  ## Title 
  $grps{Xbase}->line(id=>'baseline', 
    x1 => $x_base, y1 => $y_base, 
    x2 => $x_base+$max_len{$chrn}/$xscale, y2 => $y_base, 
  ); 
  $grps{Xbase}->text(
    x => $x_base, y => $y_base+40, 
    -cdata => "$chrn ( Drawn for $covF )", 
    'text-anchor' => "start", 
  );
  
  
  ## Base Line X. 
  my $m1 = 1000000; 
  for (my $i=0; $i*$m1<$max_len{$chrn}; $i+=1) {
    $grps{Xbase}->line(
      x1 => $x_base+$i*$m1/$xscale, y1 => $y_base, 
      x2 => $x_base+$i*$m1/$xscale, y2 => $y_base+5, 
    ); 
    $grps{Xbase}->text (
      x => $x_base+$i*$m1/$xscale, y => $y_base+20, 
      -cdata => "${i}M", 
    ); 
  }
  ## SNP sites lines. 
  my $maxY = 0; 
  for (my $i=0; $i<@{$hash{$chrn}}; $i++) {
    my ($px0, $py0, $px1, $py1); 
    # for my $grpid (qw/homoPI hete homo97/) {
    for my $grpid (qw/hete homoPI homo97/) { # 20121127 Edit the order to plot. 
# $grpid eq 'homoPI' or next; 
      if (defined $px0) {
	$py0 = $py1; 
	$py1 += $hash{$chrn}[$i][$grpCol{$grpid}]; 
      }else{
        $px0 = $px1 = $hash{$chrn}[$i][$grpCol{WinS}]; 
        $py0 = 0; 
        $py1 = $hash{$chrn}[$i][$grpCol{$grpid}]; 
      }
      $grps{$grpid}->rectangle(
        x => $x_base+$px0/$xscale, y => $y_base-$py1/$yscale, 
        width => $colWid , height=> $hash{$chrn}[$i][$grpCol{$grpid}]/$yscale, 
      ); 
    }
    $maxY < $py1 and $maxY = $py1; 
  }
  # $chrn eq 'Chr2' and die "maxY=$maxY"; 
  ## For Base Line Y 
  $grps{Ybase}->line(
    x1 => $x_base-5, y1 => $y_base, 
    x2 => $x_base-5, y2 => $y_base-$maxY/$yscale, 
  ); 
  for (my $i=0; $i <= $maxY/$yscale; $i += 50) {
    my $num = int($i*$yscale); 
    my $ty = $num/$yscale; 
    $grps{Ybase}->line(
      x1 => $x_base-5, y1 => $y_base-$ty, 
      x2 => $x_base-10, y2 => $y_base-$ty, 
    ); 
    $grps{Ybase}->text(
      x => $x_base-13, y => $y_base-$ty, 
      -cdata => "$num", 
    ); 
  }
  
  open OS,'>',"${tag}_$chrn.svg" or die; 
  print OS $svg->xmlify; 
  close OS; 
  print STDERR "${tag}_$chrn.svg\n"; 
}














