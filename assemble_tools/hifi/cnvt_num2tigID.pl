#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and die "perl $0 example_tigID in_num_file > out_tigID_file\n"; 

my $sample = shift; 

$sample =~ m!^(tig\d+)$! or die "bad example tigID [$sample]\n"; 

my $nLen = length($sample) - 3; 

while (<>) {
  chomp; 
  my @ta=split(/\t/, $_); 
  if (m!^\s*(#|$)!) {
    print "$ta[0]\n"; 
    next; 
  }
  $ta[0] =~ m!^\d+$! or die "bad number [$ta[0]]\n"; 
  my $newID = sprintf("%0${nLen}d", $ta[0]); 
  $newID = "tig$newID"; 
  print "$newID\n"; 
}

