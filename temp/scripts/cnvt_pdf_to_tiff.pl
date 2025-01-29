#!/usr/bin/perl
use strict;
use warnings;
use LogInforSunhh;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
  "out_format:s", # tiff.
  "out_dpi:i",  # 300;
  "help!"
);

$opts{'out_format'} //= 'tiff';
$opts{'out_dpi'}  //= 300;

my $help_txt = <<HH;
######################################################################
perl $0  input.pdf

-out_format       [$opts{'out_format'}] Can also be 'tiff,png'
-out_dpi          [$opts{'out_dpi'}]
######################################################################
HH

!@ARGV and &LogInforSunhh::usage($help_txt);
defined $opts{'help'} and &LogInforSunhh::usage($help_txt);

my $f=shift;
$f =~ m!^(\S+)\.(pdf|svg)$!i or die "|$f|\n";
my $baseName = $1;

for my $ofmt (split(/,/, $opts{'out_format'})) {
  $ofmt =~ s!\s!!g;
  &runCmd("convert -density $opts{'out_dpi'} $f -quality 100 $1.$ofmt");
}
# &runCmd("pdf2svg $f $1.svg");
# pdf2svg domestication_time-CLC_CLV.csv_splittime.histo.pdf domestication_time-CLC_CLV.csv_splittime.histo.svg

