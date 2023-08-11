#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
  chomp;
  m!^#! and do { print "$_\n"; next; };
  my @ta=split(/\t/, $_);
  $ta[3] =~ m![nN]! and next;
  $ta[4] ne '<INV>' and $ta[4] =~ m![nN]! and next;
  print "$_\n";
}
