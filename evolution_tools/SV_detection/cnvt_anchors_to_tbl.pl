#!/usr/bin/perl
use strict;
use warnings;

my %blk;
while (<>) {
  chomp;
  m!^\s*(#|$)! and next;
  my @ta=split(/\t/, $_);
  $ta[1] eq 'referenceStart' and next;
  $blk{$ta[8]}{'rID'} //= $ta[0];
  $blk{$ta[8]}{'qID'} //= $ta[3];
  $blk{$ta[8]}{'str'} //= $ta[6];
  $blk{$ta[8]}{'Rstart'} //= $ta[1]; $blk{$ta[8]}{'Rstart'} > $ta[1] and $blk{$ta[8]}{'Rstart'} = $ta[1];
  $blk{$ta[8]}{'Rend'}   //= $ta[2]; $blk{$ta[8]}{'Rend'}   < $ta[2] and $blk{$ta[8]}{'Rend'}   = $ta[2];
  $blk{$ta[8]}{'Qstart'} //= $ta[4]; $blk{$ta[8]}{'Qstart'} > $ta[4] and $blk{$ta[8]}{'Qstart'} = $ta[4];
  $blk{$ta[8]}{'Qend'}   //= $ta[5]; $blk{$ta[8]}{'Qend'}   < $ta[5] and $blk{$ta[8]}{'Qend'}   = $ta[5];
}
for (sort {$blk{$a}{'rID'} cmp $blk{$b}{'rID'} || $blk{$a}{'Rstart'} <=> $blk{$b}{'Rstart'} } keys %blk) {
  $blk{$_}{'lenQ'} = $blk{$_}{'Qend'} - $blk{$_}{'Qstart'} + 1;
  $blk{$_}{'lenR'} = $blk{$_}{'Rend'} - $blk{$_}{'Rstart'} + 1;
  print STDOUT join("\t", @{$blk{$_}}{qw/qID Qstart Qend str rID Rstart Rend lenQ lenR/}, , $_)."\n";
}
