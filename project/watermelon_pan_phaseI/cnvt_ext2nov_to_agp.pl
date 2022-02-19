#!/usr/bin/perl
use strict;
use warnings;

!@ARGV and die "perl $0 CM.ext2nov.tbl > CM.ext2nov.agp\n";

while (<>) {
  chomp;
  my @ta=split(/\t/, $_);
  print join("\t", $ta[4], 1, $ta[2]-$ta[1]+1, 1, "W", $ta[0], $ta[1], $ta[2], $ta[3])."\n";
}

#==> CM.ori2ext.tbl <==
#NODE_1015	1	47571	+	NODE_1015__1_47571_ext	WM1147_PI164248
#NODE_10412	1	10272	+	NODE_10412__1_10272_ext	WM1147_PI164248
#NODE_10828	1	9780	+	NODE_10828__1_9780_ext	WM1147_PI164248
#
#==> CM.ext2nov.tbl <==
#NODE_1015__1_47571_ext	16875	17446	+	WM1147_PI164248_NODE_1015_16875-17446	WM1147_PI164248
#NODE_10412__1_10272_ext	5450	6631	+	WM1147_PI164248_NODE_10412_5450-6631	WM1147_PI164248
#NODE_10828__1_9780_ext	352	8977	+	WM1147_PI164248_NODE_10828_352-8977	WM1147_PI164248
#     AGP : WM97pbV1_Chr06  1       29507460        1       W       ClaScf_0005     1       29507460        -       Scaffold5

