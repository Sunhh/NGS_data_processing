#!/usr/bin/perl
use strict;
use warnings;

my $f=shift;
open F,'<',"$f" or die;
my $is=0;
while (<F>) {
  m!^@! and next;
  chomp;
  my @ta=split(/\t/, $_);
  $ta[5] =~ m!^(\d+H)?\d+[ID]! and do {$is = 1; last;};
  $ta[5] =~ m!\d+[ID](\d+H)?$! and do {$is = 1; last;};
}
close F;
print join("\t", $f, $is)."\n";

