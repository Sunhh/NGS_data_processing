#!/usr/bin/perl
# [5/10/2022] Replace unicode characters to space. At first I need to provide a list of good characters by "perl -e 'while (<>) { split(//, $_);}'"
use strict;
use warnings;

!@ARGV and die "perl $0 good_char.list input.file > input.replaced\n";

my $f1 = shift;
open F1,'<',"$f1" or die;
my %h;
while (<F1>) {
  chomp;
  $_ = "AA${_}AA";
  my @ta=split(/\t/, $_);
  $ta[0] =~ s!^AA!!;
  $ta[-1] =~ s!AA$!!;
  $h{$ta[0]} = $ta[0];
}
close F1;
$h{"\t"} = "\t";

while (<>) {
  chomp;
  my @ta=split(//, $_);
  for my $tb (@ta) {
    defined $h{$tb} or $tb = " ";
  }
  print join("", @ta)."\n";
}

