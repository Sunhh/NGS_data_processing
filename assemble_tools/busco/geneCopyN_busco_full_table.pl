#!/usr/bin/perl
use strict; 
use warnings; 

!@ARGV and -t and die "perl $0 full_table.tsv > full_table.tsv.cnt\n"; 

my %h; 
while (<>) {
  m!^\s*#! and next; 
  m!^\s*$! and next; 
  chomp; 
  my @ta=split(/\t/, $_); 
  if ($ta[1] eq "Complete" or $ta[1] eq "Duplicated") {
    $h{"complete"}{$ta[0]} ++; 
  } elsif ($ta[1] eq "Fragmented") {
    $h{"fragmented"}{$ta[0]} ++; 
  } elsif ($ta[1] eq "Missing") {
    $h{"missing"}{$ta[0]} ++; 
  } else {
    die "unknown tag [$ta[1]]\n"; 
  }
}
for my $k1 (qw/complete fragmented missing/) {
  $h{$k1} //= {}; 
  for my $k2 (sort keys %{$h{$k1}}) {
    print join("\t", $k2, $k1, $h{$k1}{$k2})."\n"; 
  }
}
