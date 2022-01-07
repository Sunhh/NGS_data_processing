#!/usr/bin/perl
use strict; 
use warnings; 

-t and !@ARGV and die "perl $0 window_size window_step in.qual.fa > in.qual.fa.tab\n"; 

my $win_size = shift; 
my $win_step = shift; 

my ($k1, @seq); 
while (<>) {
  chomp;
  if (m!^\s*\>\s*(\S+)!) {
    &output(); 
    $k1 = $1; 
    @seq = (); 
  } else {
    s!^\s+!!g; 
    s!\s+$!!g; 
    m!^$! and next; 
    my @ta=split(/\s+/, $_); 
    push(@seq, @ta); 
  }
}

&output(); 
$k1 = ""; 
@seq = (); 

sub output {
  scalar(@seq) == 0 and return;
  $k1 eq "" and return;
  my @outArr; 
  my $i=-1; 
  for my $vdep (@seq) {
    $i++; 
    my $p=$i+1; 
    my ($i_s, $i_e); # i is number of window steps. 
    if ($p <= $win_size) {
      $i_s = 0; 
    } else {
      $i_s = int(($p-$win_size-1)/$win_step)+1; 
    }
    $i_e = int(($p-1)/$win_step); 
    for (my $j=$i_s; $j<=$i_e; $j++) {
      $outArr[$j][0] += $vdep; # total sum; 
      $outArr[$j][1] ++; # sample size: Number of numbers.
      $outArr[$j][2] //= $vdep; # minimum.
      $outArr[$j][2] > $vdep and $outArr[$j][2] = $vdep;
      $outArr[$j][3] //= $vdep; # maximum.
      $outArr[$j][3] < $vdep and $outArr[$j][3] = $vdep;
    }
  }
  $i=-1; 
  for my $o1 (@outArr) {
    $i++; 
    $o1->[4] = -1; 
    $o1->[1] > 0 and $o1->[4] = sprintf("%0.1f", $o1->[0]/$o1->[1]);
    my $p_s = $i * $win_step + 1;
    my $p_e = $p_s + $win_size -1; 
    print join("\t", $k1, $p_s, $p_e, $o1->[4], $o1->[2], $o1->[3], $o1->[1], $o1->[0])."\n"; 
  }
}

