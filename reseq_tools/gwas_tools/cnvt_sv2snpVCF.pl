#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
  m!^#! and do {print; next;};
  chomp;
  my @ta=split(/\t/, $_);
  $ta[3] = "A";
  $ta[4] =~ m!,! and die "[Err] multi-allelic input in line: $_\n";
  $ta[4] = "T";
  $ta[7] = ".";
  print join("\t", @ta)."\n";
}
