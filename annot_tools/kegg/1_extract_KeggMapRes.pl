#!/usr/bin/perl
# Copy the returned text from kegg mapper - search pathway (https://www.kegg.jp/kegg/tool/map_pathway1.html)
# Run this script to get a table with repeated KO IDs. 
use strict; 
use warnings; 

my %h; 
while (<>) {
  chomp; 
  m!^\s*$! and next; 
  if (m!^map\d+!) {
    %h = (); 
    m!^(map\d+) (.+) \(\d+\)$! or die "$_\n"; 
    $h{'mapID'} = $1; 
    $h{'mapDesc'} = $2; 
  } elsif (m!^\s\sko:(K\d+) (.+)$!) {
    my ($kid, $kdesc) = ($1, $2); 
    print STDOUT join("\t", $kid, $kdesc, $h{'mapID'}, $h{'mapDesc'})."\n";
  } else {
    die "$_\n"; 
  }
}
