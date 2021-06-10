#!/usr/bin/perl
use strict; 
use warnings; 

my %h; 
while (<>) {
  chomp; 
  my @ta=split(/\t/, $_); 
  $h{$ta[0]}{'desc'}{$ta[1]}++; 
  $h{$ta[0]}{'map'}{"$ta[2]__($ta[3])"} ++; 
  $h{$ta[0]}{'rank'} //= $.; 
}

for my $id (sort keys %h) {
  scalar(keys %{$h{$id}{'desc'}}) == 1 or die "$id\n"; 
  my ($desc) = (keys %{$h{$id}{'desc'}});
  print STDOUT join("\t", $id, $desc, join(";;", sort keys %{$h{$id}{'map'}}))."\n"; 
}
